import numpy as np
import gurobipy as gp
from gurobipy import GRB
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os

class MatrixSeriation:
    def __init__(self, matrix, file="test", method='general', time_limit=3600, extra_constraints=None, symmetric_ordering=True, output=0, eps_neigh=1.0):
        """
        Initialize the MatrixSeriation class.
        
        Parameters:
        - matrix: np.array, the matrix to be seriated
        - method: str, the optimization method to use ('linear_ordering', 'quadratic_assignment', etc.)
        - time_limit: int, the maximum allowed time for the solver in seconds
        - extra_constraints: function, a function that adds extra constraints to the model
        - enforce_square: bool, whether to enforce the matrix to be square for certain models
        - enforce_symmetric_ordering: bool, whether to enforce row and column sorting in the same order
        """
        self.matrix = np.array(matrix)
        self.file=file
        self.method = method
        self.time_limit = time_limit
        self.extra_constraints = extra_constraints
        self.symmetric_ordering = symmetric_ordering
        self.nrows, self.ncols = self.matrix.shape  # Rows and columns
        if self.nrows == self.ncols:
            self.square=True
        self.solution_info = {}
        self.output=output
        self.eps_neigh=eps_neigh
    
    def lproot(self, model, where):
        """
        Gurobi callback function to add cutting planes at MIPNODE and MIPSOL nodes.
        """
        if where == GRB.Callback.MIPNODE:  # Node solution
            nodecnt = model.cbGet(GRB.Callback.MIPNODE_NODCNT)
            if nodecnt==0:
                model._lp=model.cbGet(GRB.Callback.MIPNODE_OBJBND)
            
    def solve(self):
        """
        Solve the seriation problem using the selected method.
        """
        if self.method == 'general':
            return self._solve_general()
        elif self.method == 's':
            return self._solve_s()
        elif self.method == 'tsp':
            if self.eps_neigh<1.4:
                return self._solve_tsp()
            else:
                return self._solve_tsp8()
        elif self.method == 'spp':
            if self.eps_neigh<1.4:
                return self._solve_spp()
            else:
                return self._solve_spp8()
        else:
            raise ValueError(f"Unknown method: {self.method}")
    
    
    
    def _solve_general(self):
        """
        Solve the seriation problem using a linear ordering formulation.
        """
        model = gp.Model("four_indices_sym")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        model.setParam("MIPFocus", 2)
        

        N= range(self.nrows)
        M= range(self.ncols)

        ###x_ij/y_ij = 1 if row/column i is sorted in k position:

        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        if self.symmetric_ordering:
            y=x

        v={}
        KK={}
        for k in N:
            for l in M:
                KK[k,l]=self.get_neighbors(k,l)
                for a in KK[k,l]:
                    v[k, l, a[0], a[1]]=model.addVar(name="v[%d,%d,%d,%d]"%(k,l,a[0], a[1]))

        obj = gp.quicksum(v[k, l, a[0], a[1]] for k in N for l in M for a in KK[k,l])
        model.setObjective(obj, GRB.MINIMIZE)
        # for i in N:
        #     for j in N:
        #         for k in N:
        #             for l in N:
        #                 model.addConstr(y[i,j,k,l] >= x[i,k]+x[j,l]-1)
        #                 model.addConstr(y[i,j,k,l] <= x[i,k])
        #                 model.addConstr(y[i,j,k,l] <= x[j,l])
        # for k in N:
        #     for l in N:
        #         for a in KK[k,l]:
        #             model.addConstr(v[k,l, a[0], a[1]] >=  gp.quicksum(self.matrix[i,j]*y[i,j,k,l] for i in N for j in N) - gp.quicksum(self.matrix[i,j]*y[i,j,a[0], a[1]] for i in N for j in N))
        #             model.addConstr(v[k,l, a[0], a[1]] >= -gp.quicksum(self.matrix[i,j]*y[i,j,k,l] for i in N for j in N) + gp.quicksum(self.matrix[i,j]*y[i,j,a[0], a[1]] for i in N for j in N))
        for k in N:
            for l in M:
                for a in KK[k,l]:
                    model.addConstr(v[k,l, a[0], a[1]] >=  gp.quicksum(self.matrix[i,j]*x[i,k]*y[j,l] for i in N for j in M) - gp.quicksum(self.matrix[i,j]*x[i,a[0]]*y[j,a[1]] for i in N for j in M))
                    model.addConstr(v[k,l, a[0], a[1]] >=  -gp.quicksum(self.matrix[i,j]*x[i,k]*y[j,l] for i in N for j in M) + gp.quicksum(self.matrix[i,j]*x[i,a[0]]*y[j,a[1]] for i in N for j in M))

        for i in N:
            model.addConstr(gp.quicksum(x[i,k] for k in N)==1)
            model.addConstr(gp.quicksum(x[k,i] for k in N)==1)



        if not self.symmetric_ordering:
            for i in M:
                model.addConstr(gp.quicksum(y[i,k] for k in M)==1)
                model.addConstr(gp.quicksum(y[k,i] for k in M)==1)

        model._lproot=0
        model.optimize(self.lproot)
        if model.Status == GRB.INFEASIBLE:
            model.computeIIS()
            model.write('iismodel_general.ilp')
            print('infeasible model')

        transformed_matrix = np.empty((self.nrows,self.ncols))
        for k in N:
            for l in N:
                transformed_matrix[k,l]=sum(self.matrix[i,j]*x[i,k].x*y[j,l].x for i in N for j in M)
    
        self.transformed_matrix=transformed_matrix
        self.order_rows=[int(sum(k*x[i,k].x for k in N)) for i in N]
        self.order_cols=[int(sum(k*y[i,k].x for k in M)) for i in M]

        self._extract_solution_info(model)

        return transformed_matrix
    
    
    def _solve_s(self):

        # Model
        model = gp.Model("four_indices_s_nonsym")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        model.setParam("MIPFocus", 2)


        N=range(self.nrows)
        M=range(self.ncols)
        # Decision variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        if self.symmetric_ordering:
            y=x

        s = model.addVars(self.nrows, self.ncols, ub=1, name="s")
        v = {}  # Dynamic variables

        # Constraints
        # Each row and column of x sums to 1
        model.addConstrs((gp.quicksum(x[i, k] for k in N) == 1 for i in N), name="RowSumX")
        model.addConstrs((gp.quicksum(x[k, i] for k in N) == 1 for i in N), name="ColSumX")

        if not self.symmetric_ordering:
            model.addConstrs((gp.quicksum(y[i, k] for k in M) == 1 for i in M), name="RowSumY")
            model.addConstrs((gp.quicksum(y[k, i] for k in M) == 1 for i in M), name="ColSumY")

        # # Order constraint
        model.addConstr(
            gp.quicksum(2**(k+1) * x[0, k] for k in N) >= gp.quicksum(2**(k+1) * x[1, k] for k in N),
            name="OrderConstraint"
        )

        # Dynamic variables and constraints
        KK={}
        for k in N:
            for l in M:
                KK[k,l]=self.get_neighbors(k,l)
                for a in KK[k,l]:
                    v[k, l, a[0], a[1]] = model.addVar(vtype=GRB.CONTINUOUS, name=f"v[{k},{l},{a[0]},{a[1]}]")                

        # s constraints
        model.addConstrs(
            s[k, l] >= x[i, k] + gp.quicksum(self.matrix[i,j]*y[j, l] for j in M if self.matrix[i,j]>0.001) - 1
            for i in N for k in N for l in M
        )
        model.addConstrs(
            s[k, l] >= y[j,l] + gp.quicksum(self.matrix[i,j]*x[i,k] for i in N if self.matrix[i,j]>0.001) - 1
            for j in M for k in N for l in M
        )
        model.addConstrs(
            s[k, l] <= 1 - x[i,k] + gp.quicksum(self.matrix[i,j]*y[j,l] for j in M if self.matrix[i,j]>0.001)
            for i in N for k in N for l in M
        )
        model.addConstrs(
            s[k, l] <= 1 - y[j,l] + gp.quicksum(self.matrix[i,j]*x[i,k] for i in N if self.matrix[i,j]>0.001)
            for j in M for k in N for l in M
        )

        # v constraints
        model.addConstrs(
            v[k, l, kp, lp] >= s[k, l] - s[kp, lp] for k, l, kp, lp in v.keys()
        )
        model.addConstrs(
            v[k, l, kp, lp] >= s[kp, lp] - s[k, l] for k, l, kp, lp in v.keys()
        )

        model.addConstr(gp.quicksum(s[k,l] for k in N for l in M) == sum(self.matrix[i,j] for i in N for j in M))

        model.addConstrs(
            gp.quicksum(s[k,l] for l in M) == gp.quicksum(self.matrix[i,j]*x[i,k] for i in N for j in M)
              for k in N
        )
        if not self.symmetric_ordering:
            model.addConstrs(
                gp.quicksum(s[k,l] for k in N) == gp.quicksum(self.matrix[i,j]*y[j,l] for i in N for j in M)
                for l in M
            )

        # Objective: Minimize sum of v variables
        model.setObjective(
            0.5*gp.quicksum(v[k, l, kp, lp] for k, l, kp, lp in v.keys()), GRB.MINIMIZE
        )

        model._lproot=0
        model.optimize(self.lproot)
        # Extract solution information
        #self._extract_solution_info(model)
        
        # transformed matrix
        transformed_matrix = np.empty((self.nrows,self.ncols))
        for k in N:
            for l in M:
                transformed_matrix[k,l]=sum(self.matrix[i,j]*x[i,k].x*y[j,l].x for i in N for j in M)
    
        self.transformed_matrix=transformed_matrix
        self.order_rows=[int(sum(k*x[i,k].x for k in N)) for i in N]
        self.order_cols=[int(sum(k*y[i,k].x for k in M)) for i in M]

        self._extract_solution_info(model)

        return transformed_matrix



    def _solve_tsp(self):

        # Modelcle
        model = gp.Model("tsp")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        model.setParam("MIPFocus", 2)


        N=range(self.nrows)
        M=range(self.ncols)
        # Compute costs
        C1={}
        C2={ }
        for i in N:
            for j in N:
                C1[i, j] =  sum(abs(self.matrix[i, k] - self.matrix[j, k]) for k in M)
        for i in M:
            for j in M:
                C2[i, j] =  sum(abs(self.matrix[k, i] - self.matrix[k, j]) for k in N)


        # Variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        

        gx = model.addVars(self.nrows, self.nrows, name="gx")
        gy = model.addVars(self.ncols, self.ncols, name="gy")

        zx = model.addVars(self.nrows, vtype=GRB.BINARY, name="zx")
        zy = model.addVars(self.ncols, vtype=GRB.BINARY, name="zy")

        if self.symmetric_ordering:
            y=x
            gy=gx
            zy=zx


        # Objective function
        model.setObjective(gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N) + gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M), GRB.MINIMIZE)


        # Constraints
        model.addConstr(zx[0] == 0)
        model.addConstr(gp.quicksum(zx[i] for i in N) == 1)
        for i in N:
            model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) + zx[i] == 1)
            x[i, i].ub=0
            model.addConstr(gp.quicksum(x[j,i] for j in N if i != j) <= 1)
            model.addConstr(gp.quicksum(gx[i, j] for j in N) >= gp.quicksum(gx[j, i] for j in N) - (self.nrows) * zx[i] + 1)
            for j in N:
                if j<i:
                    model.addConstr(x[i, j] + x[j, i] <= 1)
                model.addConstr(gx[i, j] <= (self.nrows-1) * x[i, j])
                #model.addConstr(gx[i, j] >=  x[i, j])
        
            
        if not self.symmetric_ordering:
            model.addConstr(zy[0] == 0)
            model.addConstr(gp.quicksum(zy[i] for i in M) == 1)
            for i in M:
                model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) + zy[i] == 1)
                model.addConstr(gp.quicksum(y[j,i] for j in M if i != j) <= 1)
                model.addConstr(y[i, i] == 0)
                model.addConstr(gp.quicksum(gy[i, j] for j in M) >= gp.quicksum(gy[j, i] for j in M) - (self.ncols) * zy[i] + 1)
                for j in M:
                    if j<i:
                        model.addConstr(y[i, j] + y[j, i] <= 1)
                    model.addConstr(gy[i, j] <= (self.ncols-1) * y[i, j])
                    #model.addConstr(gy[i, j] >=  y[i, j])


        # Solve model
        model._lproot=0
        model.optimize(self.lproot)
        #model.write("tsp.sol")

        # Extract solution
        self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M)-1 if zy[i].x < 0.5 else self.ncols-1 for i in M])
        self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N)-1 if zx[i].x < 0.5 else self.nrows-1 for i in N])
        
        # for k in range(1,self.nrows):
        #     for i in N:
        #         for j in N:
        #             if gx[i,j].x==k:
        #                 print(f"{k}: {i}->{j}")
        # print(self.order_rows, self.order_cols)

        transformed_matrix = self.matrix[self.order_rows][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix

    def _solve_tsp8(self):

        # Modelcle
        model = gp.Model("tsp")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        model.setParam("MIPFocus", 2)


        N=range(self.nrows)
        M=range(self.ncols)
        # Compute costs
        C1={}
        C2={ }
        for i in N:
            for j in N:
                C1[i, j] =  sum(abs(self.matrix[i, k] - self.matrix[j, k]) for k in M)
        for i in M:
            for j in M:
                C2[i, j] =  sum(abs(self.matrix[k, i] - self.matrix[k, j]) for k in N)
        C={}
        for i1 in N:
            for i2 in N:
                for j1 in M:
                    for j2 in M:
                        C[i1,i2,j1,j2]= np.abs(self.matrix[i1,j1]-self.matrix[i2,j2])+\
                        np.abs(self.matrix[i1,j2]-self.matrix[i2,j1])


        # Variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        

        gx = model.addVars(self.nrows, self.nrows, name="gx")
        gy = model.addVars(self.ncols, self.ncols, name="gy")

        zx = model.addVars(self.nrows, vtype=GRB.BINARY, name="zx")
        zy = model.addVars(self.ncols, vtype=GRB.BINARY, name="zy")

        if self.symmetric_ordering:
            y=x
            gy=gx
            zy=zx

        C={}
        for i1 in N:
            for i2 in N:
                for j1 in M:
                    for j2 in M:
                        C[i1,i2,j1,j2]= np.abs(self.matrix[i1,j1]-self.matrix[i2,j2])+\
                        np.abs(self.matrix[i1,j2]-self.matrix[i2,j1])
                        

                        
       # Objective function
        obj1=gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N)
        obj2=gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M)
        obj12=gp.quicksum(C[i1,i2,j1,j2]*x[i1,i2]*y[j1,j2] for i1 in N for i2 in N for j1 in M for j2 in M)
        model.setObjective(obj1+obj2+obj12, GRB.MINIMIZE)


        # Constraints
        model.addConstr(zx[0] == 0)
        model.addConstr(gp.quicksum(zx[i] for i in N) == 1)
        for i in N:
            model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) + zx[i] == 1)
            x[i, i].ub=0
            model.addConstr(gp.quicksum(x[j,i] for j in N if i != j) <= 1)
            model.addConstr(gp.quicksum(gx[i, j] for j in N) >= gp.quicksum(gx[j, i] for j in N) - (self.nrows) * zx[i] + 1)
            for j in N:
                if j<i:
                    model.addConstr(x[i, j] + x[j, i] <= 1)
                model.addConstr(gx[i, j] <= (self.nrows-1) * x[i, j])
                #model.addConstr(gx[i, j] >=  x[i, j])
        
            
        if not self.symmetric_ordering:
            model.addConstr(zy[0] == 0)
            model.addConstr(gp.quicksum(zy[i] for i in M) == 1)
            for i in M:
                model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) + zy[i] == 1)
                model.addConstr(gp.quicksum(y[j,i] for j in M if i != j) <= 1)
                model.addConstr(y[i, i] == 0)
                model.addConstr(gp.quicksum(gy[i, j] for j in M) >= gp.quicksum(gy[j, i] for j in M) - (self.ncols) * zy[i] + 1)
                for j in M:
                    if j<i:
                        model.addConstr(y[i, j] + y[j, i] <= 1)
                    model.addConstr(gy[i, j] <= (self.ncols-1) * y[i, j])
                    #model.addConstr(gy[i, j] >=  y[i, j])


        # Solve model
        model._lproot=0
        model.optimize(self.lproot)
        #model.write("tsp.sol")

        # Extract solution
        self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M)-1 if zy[i].x < 0.5 else self.ncols-1 for i in M])
        self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N)-1 if zx[i].x < 0.5 else self.nrows-1 for i in N])
        
        # for k in range(1,self.nrows):
        #     for i in N:
        #         for j in N:
        #             if gx[i,j].x==k:
        #                 print(f"{k}: {i}->{j}")
        # print(self.order_rows, self.order_cols)

        transformed_matrix = self.matrix[self.order_rows][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix

    def _solve_spp(self):

        # Modelcle
        model = gp.Model("tsp")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        model.setParam("MIPFocus", 2)


        N=range(self.nrows)
        M=range(self.ncols)
        # Compute costs
        C1={}
        C2={ }
        for i in N:
            for j in N:
                C1[i, j] =  sum(abs(self.matrix[i, k] - self.matrix[j, k]) for k in M)
        for i in M:
            for j in M:
                C2[i, j] =  sum(abs(self.matrix[k, i] - self.matrix[k, j]) for k in N)


        # Variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        
        ux = model.addVars(self.nrows, lb=1, ub=self.nrows, name="ux")
        uy = model.addVars(self.ncols,lb=1, ub=self.ncols, name="sy")

        sx = model.addVars(self.nrows, vtype=GRB.BINARY,name="sx")
        sy = model.addVars(self.ncols, vtype=GRB.BINARY,name="sy")
        tx = model.addVars(self.nrows, vtype=GRB.BINARY,name="tx")
        ty = model.addVars(self.ncols, vtype=GRB.BINARY,name="ty")


        if self.symmetric_ordering:
            y=x
            sy=sx
            ty=tx
            uy=ux


        # Objective function
        model.setObjective(gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N) + gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M), GRB.MINIMIZE)


        # Constraints
        model.addConstr(sx[0] == 0)
        model.addConstr(gp.quicksum(sx[i] for i in N) == 1)
        model.addConstr(gp.quicksum(tx[i] for i in N) == 1)
        for i in N:
            model.addConstr(sx[i]+tx[i] <=1)
            model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) == 1-tx[i])
            x[i, i].ub=0
            model.addConstr(gp.quicksum(x[j,i] for j in N if i != j) == 1-sx[i])
            for j in N:
                model.addConstr(ux[i]-ux[j] + self.nrows*x[i,j] <= self.nrows-1)
                if j<i:
                    model.addConstr(x[i, j] + x[j, i] <= 1)
        
            
        if not self.symmetric_ordering:
            model.addConstr(sy[0] == 0)
            model.addConstr(gp.quicksum(sy[i] for i in M) == 1)
            model.addConstr(gp.quicksum(ty[i] for i in M) == 1)
            for i in M:
                model.addConstr(sy[i]+ty[i] <=1)
                model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) == 1-ty[i])
                model.addConstr(gp.quicksum(y[j, i] for j in M if i != j)  == 1-sy[i])
                model.addConstr(y[i, i] == 0)
                for j in M:
                    model.addConstr(uy[i]-uy[j] + self.ncols*y[i,j] <= self.ncols-1)
                    if j<i:
                        model.addConstr(y[i, j] + y[j, i] <= 1)


        # Solve model
        model._lproot=0
        model.optimize(self.lproot)
        #model.write("spp.sol")

        # Extract solution
        orderx=[]
        for i in N:
            if sx[i].x>0.5:
                orderx.append(i)
                i0=i
        while tx[i0].x<0.5:
            for j in N:
                if x[i0,j].x>0.5:
                    orderx.append(j)
                    i0=j
                    break

        if self.symmetric_ordering:
            ordery=orderx
        else:
            ordery=[]
            for i in M:
                if sy[i].x>0.5:
                    ordery.append(i)
                    i0=i
            while ty[i0].x<0.5:
                for j in M:
                    if y[i0,j].x>0.5:
                        ordery.append(j)
                        i0=j
                        break




        self.order_cols = ordery
        self.order_rows = orderx
        
        # for k in range(1,self.nrows):
        #     for i in N:
        #         for j in N:
        #             if gx[i,j].x==k:
        #                 print(f"{k}: {i}->{j}")
        # print(self.order_rows, self.order_cols)

        transformed_matrix = self.matrix[self.order_rows][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix
    
    def _solve_spp8(self):

        # Modelcle
        model = gp.Model("tsp")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        model.setParam("MIPFocus", 2)


        N=range(self.nrows)
        M=range(self.ncols)
        # Compute costs
        C1={}
        C2={ }
        for i in N:
            for j in N:
                C1[i, j] =  sum(abs(self.matrix[i, k] - self.matrix[j, k]) for k in M)
        for i in M:
            for j in M:
                C2[i, j] =  sum(abs(self.matrix[k, i] - self.matrix[k, j]) for k in N)


        # Variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        
        ux = model.addVars(self.nrows, lb=1, ub=self.nrows, name="ux")
        uy = model.addVars(self.ncols,lb=1, ub=self.ncols, name="sy")

        sx = model.addVars(self.nrows, vtype=GRB.BINARY,name="sx")
        sy = model.addVars(self.ncols, vtype=GRB.BINARY,name="sy")
        tx = model.addVars(self.nrows, vtype=GRB.BINARY,name="tx")
        ty = model.addVars(self.ncols, vtype=GRB.BINARY,name="ty")


        if self.symmetric_ordering:
            y=x
            sy=sx
            ty=tx
            uy=ux


        # Objective function
        C={}
        for i1 in N:
            for i2 in N:
                for j1 in M:
                    for j2 in M:
                        C[i1,i2,j1,j2]= np.abs(self.matrix[i1,j1]-self.matrix[i2,j2])+\
                        np.abs(self.matrix[i1,j2]-self.matrix[i2,j1])
                        

                        
       # Objective function
        obj1=gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N)
        obj2=gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M)
        obj12=gp.quicksum(C[i1,i2,j1,j2]*x[i1,i2]*y[j1,j2] for i1 in N for i2 in N for j1 in M for j2 in M)
        model.setObjective(obj1+obj2+obj12, GRB.MINIMIZE)


        # Constraints
        model.addConstr(sx[0] == 0)
        model.addConstr(gp.quicksum(sx[i] for i in N) == 1)
        model.addConstr(gp.quicksum(tx[i] for i in N) == 1)
        for i in N:
            model.addConstr(sx[i]+tx[i] <=1)
            model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) == 1-tx[i])
            x[i, i].ub=0
            model.addConstr(gp.quicksum(x[j,i] for j in N if i != j) == 1-sx[i])
            for j in N:
                model.addConstr(ux[i]-ux[j] + self.nrows*x[i,j] <= self.nrows-1)
                if j<i:
                    model.addConstr(x[i, j] + x[j, i] <= 1)
        
            
        if not self.symmetric_ordering:
            model.addConstr(sy[0] == 0)
            model.addConstr(gp.quicksum(sy[i] for i in M) == 1)
            model.addConstr(gp.quicksum(ty[i] for i in M) == 1)
            for i in M:
                model.addConstr(sy[i]+ty[i] <=1)
                model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) == 1-ty[i])
                model.addConstr(gp.quicksum(y[j, i] for j in M if i != j)  == 1-sy[i])
                model.addConstr(y[i, i] == 0)
                for j in M:
                    model.addConstr(uy[i]-uy[j] + self.ncols*y[i,j] <= self.ncols-1)
                    if j<i:
                        model.addConstr(y[i, j] + y[j, i] <= 1)


        # Solve model
        model._lproot=0
        model.optimize(self.lproot)
        #model.write("spp.sol")

        # Extract solution
        orderx=[]
        for i in N:
            if sx[i].x>0.5:
                orderx.append(i)
                i0=i
        while tx[i0].x<0.5:
            for j in N:
                if x[i0,j].x>0.5:
                    orderx.append(j)
                    i0=j
                    break

        if self.symmetric_ordering:
            ordery=orderx
        else:
            ordery=[]
            for i in M:
                if sy[i].x>0.5:
                    ordery.append(i)
                    i0=i
            while ty[i0].x<0.5:
                for j in M:
                    if y[i0,j].x>0.5:
                        ordery.append(j)
                        i0=j
                        break




        self.order_cols = ordery
        self.order_rows = orderx
        
        # for k in range(1,self.nrows):
        #     for i in N:
        #         for j in N:
        #             if gx[i,j].x==k:
        #                 print(f"{k}: {i}->{j}")
        # print(self.order_rows, self.order_cols)

        transformed_matrix = self.matrix[self.order_rows][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix
    
    def _compute_error_matrix(self, A):
        
        error=np.zeros((self.nrows, self.ncols))
        for i in range(self.nrows):
            for j in range(self.ncols):
                error[i,j]=sum(np.abs(A[i,j]-A[a[0],a[1]]) for a in self.get_neighbors(i,j))

            
        return 0.5*sum(sum(error))



    def _extract_solution_info(self, model):
        """Extracts relevant solution information after solving the model."""
        self.solution_info = {
            "file": self.file,
            "Formulation": self.method,
            "nrows": self.nrows,
            "ncols": self.ncols,
            "Sym": self.symmetric_ordering,
            "eps_neigh": self.eps_neigh,
            "CPUTime": model.Runtime,
            "ObjectiveValue": model.ObjVal if model.status != GRB.INFEASIBLE else None,
            "MIPGap": model.MIPGap,
            "LPRelaxation at Root": model._lproot,
            "Work": model.getAttr(GRB.Attr.Work),
            "NodeCount": model.NodeCount,
            "Memory": model.getAttr(GRB.Attr.MemUsed),
            "NumBinVars": model.NumBinVars,
            "NumCtrs": model.NumConstrs,
            "Homogeneity Original": self.Homogeneity(self.matrix),
            "Homogeneity Transformed": self.Homogeneity(self.transformed_matrix),
            "Order Rows": self.order_rows,
            "Order Cols": self.order_cols,
            "Error Original": self._compute_error_matrix(self.matrix),
            "Error Transformed": self._compute_error_matrix(self.transformed_matrix)
            #"TMatrix":  self.transformed_matrix
        }

    def get_neighbors(self, i, j):
        """Computes neighbors for a given (i, j) in the matrix."""
        neighbors=[]
        for i1 in range(self.nrows):
            for j1 in range(self.ncols):
                if 0.1 <= np.linalg.norm([i-i1, j-j1])<=self.eps_neigh:
                    neighbors.append((i1,j1))

        return neighbors
    
    def get_adjacent_indices(self, i, j):
        """Computes adjacent indices for a given (i, j) in the matrix."""
        adjacent_indices = []
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:  # Up, Down, Left, Right
            ni, nj = i + di, j + dj
            if 0 <= ni < self.nrows and 0 <= nj < self.ncols:
                adjacent_indices.append((ni, nj))
        return adjacent_indices
    
    def get_adjacent_rows(self, i):
        """Computes adjacent indices for a given (i, j) in the matrix."""
        adjacent_rows = []
        for di in [-1, 1]:  # Up, Down
            ni = i + di
            if 0 <= ni < self.nrows:
                adjacent_rows.append(ni)
        return adjacent_rows
    
    def get_adjacent_columns(self, j):
        """Computes adjacent indices for a given (i, j) in the matrix."""
        adjacent_columns = []
        for dj in [-1, 1]:  # Left, Right
            nj = j + dj
            if 0 <= nj < self.ncols:
                adjacent_columns.append(nj)
        return adjacent_columns
    
    def Homogeneity(self, matrix):

        N=range(self.nrows)
        M=range(self.ncols)
        HOM=np.empty((self.nrows,self.ncols))
        for i in N:
            for j in M:
                HOM[i,j]=(1/len(self.get_neighbors(i,j)))*sum(abs(matrix[i,j]-matrix[i1,j1]) for (i1,j1) in self.get_adjacent_indices(i,j))

        return (1/(self.nrows*self.ncols))*sum(sum(HOM))


    def plot_matrices(self):
        """Plots the original and transformed matrices side by side."""
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        cmap = sns.color_palette("coolwarm", as_cmap=True)
        axisAx=[i+1 for i in range(self.ncols)]
        axisAy=[i+1 for i in range(self.nrows)]
        vmin=self.matrix.min()
        vmax=self.matrix.max()
        heat1=sns.heatmap(self.matrix, annot=False, fmt=".2f", cmap=cmap, xticklabels=axisAx, yticklabels=axisAy, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[0], annot_kws={"size": 6})
        heat1.invert_yaxis()
        hom1=self.Homogeneity(self.matrix)
        axes[0].set_title(f"Original Matrix - {hom1}")
        #axes[0].axis('equal')
        axisBx=self.order_cols
        axisBy=self.order_rows
        heat2=sns.heatmap(self.transformed_matrix, annot=False, fmt=".2f", xticklabels=axisBx, yticklabels=axisBy, cmap=cmap, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[1], annot_kws={"size": 6})
        heat2.invert_yaxis()
        hom2=self.Homogeneity(self.transformed_matrix)
        axes[1].set_title(f"Transformed Matrix - {hom2}")
        #axes[1].axis('equal')

        plt.savefig(f"{self.file[:-4]}_{self.method}_sym{self.symmetric_ordering}_eps{self.eps_neigh}.png", dpi=300, bbox_inches='tight')
        plt.close()
    
def generate_mosel_style_matrix(n, m):
    """
    Generate a structured matrix based on the Mosel code, apply random swaps.
    
    Parameters:
    - n (int): Number of rows.
    - m (int): Number of columns.
    
    Returns:
    - a (ndarray): The generated matrix.
    """
    
    # Step 1: Create an empty matrix
    a = np.zeros((n, m), dtype=int)
    
    # Step 2: Fill the matrix with the Mosel-defined pattern
    a[:n//2, :m//2] = 1  # Top-left quadrant
    a[n//2:, m//2:] = 1  # Bottom-right quadrant
    #a[n//2:, :m//2] = 3  # Bottom-right quadrant

    # Step 3: Randomly swap rows 100 times
    for _ in range(100):
        ii, iii = np.random.randint(0, n, 2)
        a[[ii, iii], :] = a[[iii, ii], :]

    # Step 4: Randomly swap columns 100 times
    for _ in range(100):
        jj, jjj = np.random.randint(0, m, 2)
        a[:, [jj, jjj]] = a[:, [jjj, jj]]

    return a

def generate_easy_instances(n: int):
    p = np.random.uniform(0, 100, n)
    return np.abs(p[:, np.newaxis] - p[np.newaxis, :])

def generate_sqr_instances(n: int):
    p = np.random.uniform(0, 100, (n,2))
    return np.sqrt(((p[:, np.newaxis, :] - p[np.newaxis, :, :]) ** 2).sum(axis=2))

def generate_nonsqr_instances(n: int, m:int):
    p = np.random.uniform(0, 100, (n,2))
    q = np.random.uniform(0, 100, (m,2))
    return np.sqrt(((p[:, np.newaxis, :] - q[np.newaxis, :, :]) ** 2).sum(axis=2))


def generate_binary_instances(n, m, density=0.5):
    # Total number of 1s to place
    num_ones = int(n * m * density)
    # Create a flat array with the correct number of 1s and 0s
    matrix = np.zeros(n * m, dtype=int)
    matrix[:num_ones] = 1  # Set the first `num_ones` elements to 1
    # Shuffle to distribute 1s randomly
    np.random.shuffle(matrix)
    # Reshape into an n x m matrix
    return matrix.reshape(n, m)

    



    