import numpy as np
import gurobipy as gp
from gurobipy import GRB
import seaborn as sns
import pandas as pd
import os

class MatrixSeriation:
    def __init__(self, matrix, file="test", method='general', time_limit=3600, extra_constraints=None, coordinated_ordering=True, output=0, eps_neigh="VN"):
        """
        Initialize the MatrixSeriation class.
        
        Parameters:
        - matrix: np.array, the matrix to be seriated
        - method: str, the optimization method to use ('linear_ordering', 'quadratic_assignment', etc.)
        - time_limit: int, the maximum allowed time for the solver in seconds
        - extra_constraints: function, a function that adds extra constraints to the model
        - enforce_square: bool, whether to enforce the matrix to be square for certain models
        - enforce_coordinated_ordering: bool, whether to enforce row and column sorting in the same order
        """
        self.matrix = np.array(matrix)
        self.file=file
        self.method = method
        self.time_limit = time_limit
        self.extra_constraints = extra_constraints
        self.coordinated_ordering = coordinated_ordering
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
            if self.eps_neigh=="VonNeumann":
                return self._solve_tsp_VonNeumann()
            elif self.eps_neigh=="Moore":
                return self._solve_tsp_Moore()
            elif self.eps_neigh == "Cross":
                return self._solve_tsp_Cross()
            elif self.eps_neigh == "ME":
                return self._solve_tsp_ME()
            else:
                return "Non supported neighborhood"
        # elif self.method == 'spp':
        #     if self.eps_neigh<1.4:
        #         return self._solve_spp()
        #     else:
        #         return self._solve_spp8()
        else:
            raise ValueError(f"Unknown method: {self.method}")
    
    
    
    def _solve_general(self):
        """
        Solve the seriation problem using a linear ordering formulation.
        """
        model = gp.Model("four_indices_sym")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("MIPFocus", 2)
        

        N= range(self.nrows)
        M= range(self.ncols)

        ###x_ij/y_ij = 1 if row/column i is sorted in k position:

        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])

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

        model.addConstr(gp.quicksum(2**k * x[0,k] for k in N) <= gp.quicksum(2**k * x[1,k] for k in N))
    

        if not self.coordinated_ordering:
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
        self.order_rows=[int(sum((i+1)*x[i,k].x for i in N)) for k in N]
        self.order_cols=[int(sum((i+1)*y[i,k].x for i in M)) for k in M]

        self._extract_solution_info(model)

        return transformed_matrix
    
    
    def _solve_s(self):

        # Model
        model = gp.Model("four_indices_s_nonsym")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("MIPFocus", 2)


        N=range(self.nrows)
        M=range(self.ncols)
        # Decision variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])

        s = model.addVars(self.nrows, self.ncols, ub=1, name="s")
        v = {}  # Dynamic variables

        # Constraints
        # Each row and column of x sums to 1
        model.addConstrs((gp.quicksum(x[i, k] for k in N) == 1 for i in N), name="RowSumX")
        model.addConstrs((gp.quicksum(x[k, i] for k in N) == 1 for i in N), name="ColSumX")

        if not self.coordinated_ordering:
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
        if not self.coordinated_ordering:
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
        self.order_rows=[int(sum((i+1)*x[i,k].x for i in N)) for k in N]
        self.order_cols=[int(sum((i+1)*y[i,k].x for i in M)) for k in M]

        self._extract_solution_info(model)

        return transformed_matrix



    def _solve_tsp_VonNeumann(self):

        # Modelcle
        model = gp.Model("tsp")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("MIPFocus", 2)


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

        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])
                    model.addConstr(gy[i,k]==gx[i,k])
                model.addConstr(zy[i]==zx[i])


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
        
            
        if not self.coordinated_ordering:
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



        # Extract solution
        self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M) if zy[i].x < 0.5 else self.ncols-1 for i in M])
        self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N) if zx[i].x < 0.5 else self.nrows-1 for i in N])
        
        transformed_matrix = self.matrix[self.order_rows,:][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix

    def _solve_tsp_Moore(self):

        # Modelcle
        model = gp.Model("tsp")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("MIPFocus", 2)


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

        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])
                    model.addConstr(gy[i,k]==gx[i,k])
                model.addConstr(zy[i]==zx[i])

        C={}
        for i1 in N:
            for i2 in N:
                for j1 in M:
                    for j2 in M:
                        C[i1,i2,j1,j2]= np.abs(self.matrix[i1,j1]-self.matrix[i2,j2])+ np.abs(self.matrix[i1,j2]-self.matrix[i2,j1])
                        

                        
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
        
            
        if not self.coordinated_ordering:
            model.addConstr(zy[0] == 0)
            model.addConstr(gp.quicksum(zy[i] for i in M) == 1)
            for i in M:
                model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) + zy[i] == 1)
                model.addConstr(gp.quicksum(y[j, i] for j in M if i != j) <= 1)
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
        self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M) if zy[i].x < 0.5 else self.ncols for i in M])
        self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N) if zx[i].x < 0.5 else self.nrows for i in N])
        
        # for k in range(1,self.nrows):
        #     for i in N:
        #         for j in N:
        #             if gx[i,j].x==k:
        #                 print(f"{k}: {i}->{j}")
        # print(self.order_rows, self.order_cols)

        transformed_matrix = self.matrix[self.order_rows,:][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix
    
    def _solve_tsp_Cross(self):

        # Modelcle
        model = gp.Model("tsp_cross")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("MIPFocus", 2)


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

        u = model.addVars(self.nrows, name="u")
        v = model.addVars(self.nrows, name="u")
        

        gx = model.addVars(self.nrows, self.nrows, name="gx")
        gy = model.addVars(self.ncols, self.ncols, name="gy")

        zx = model.addVars(self.nrows, vtype=GRB.BINARY, name="zx")
        zy = model.addVars(self.ncols, vtype=GRB.BINARY, name="zy")

        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])
                    model.addConstr(gy[i,k]==gx[i,k])
                model.addConstr(zy[i]==zx[i])


                        
       # Objective function
        obj1=gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N)
        obj2=gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M)
        obj12=gp.quicksum(u[i] for i in N) + gp.quicksum(v[j] for j in M)
        model.setObjective(obj1+obj2+obj12, GRB.MINIMIZE)

        for i in N:
            Delta1=np.max([C1[i,j] for j in N])
            model.addConstr(u[i] >= gp.quicksum(C1[i,l]*x[k,l] for l in N) + Delta1*(1-x[i,k]))
        for j in M:
            Delta2=np.max([C2[j,l] for l in M])
            model.addConstr(v[j] >= gp.quicksum(C2[j,l]*y[k,l] for l in M) + Delta2*(1-x[j,k]))

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
        
            
        if not self.coordinated_ordering:
            model.addConstr(zy[0] == 0)
            model.addConstr(gp.quicksum(zy[i] for i in M) == 1)
            for i in M:
                model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) + zy[i] == 1)
                model.addConstr(gp.quicksum(y[j, i] for j in M if i != j) <= 1)
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
        self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M) if zy[i].x < 0.5 else self.ncols for i in M])
        self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N) if zx[i].x < 0.5 else self.nrows for i in N])
        
        # for k in range(1,self.nrows):
        #     for i in N:
        #         for j in N:
        #             if gx[i,j].x==k:
        #                 print(f"{k}: {i}->{j}")
        # print(self.order_rows, self.order_cols)

        transformed_matrix = self.matrix[self.order_rows,:][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix
    
    def _solve_tsp_ME(self):

        # Modelcle
        model = gp.Model("tsp_cross")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("MIPFocus", 2)


        N=range(self.nrows)
        M=range(self.ncols)
        # Compute costs
        C1={}
        C2={ }
        for i in N:
            for j in N:
                C1[i, j] =  sum(self.matrix[i, k]*self.matrix[j, k] for k in M)
        for i in M:
            for j in M:
                C2[i, j] =  sum(self.matrix[k, i]*self.matrix[k, j] for k in N)


        # Variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")

        gx = model.addVars(self.nrows, self.nrows, name="gx")
        gy = model.addVars(self.ncols, self.ncols, name="gy")

        zx = model.addVars(self.nrows, vtype=GRB.BINARY, name="zx")
        zy = model.addVars(self.ncols, vtype=GRB.BINARY, name="zy")

        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])
                    model.addConstr(gy[i,k]==gx[i,k])
                model.addConstr(zy[i]==zx[i])


                        
       # Objective function
        obj1=gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N)
        obj2=gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M)
        model.setObjective(obj1+obj2, GRB.MAXIMIZE)


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
        
            
        if not self.coordinated_ordering:
            model.addConstr(zy[0] == 0)
            model.addConstr(gp.quicksum(zy[i] for i in M) == 1)
            for i in M:
                model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) + zy[i] == 1)
                model.addConstr(gp.quicksum(y[j, i] for j in M if i != j) <= 1)
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
        self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M) if zy[i].x < 0.5 else self.ncols-1 for i in M])
        self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N) if zx[i].x < 0.5 else self.nrows-1 for i in N])
        
        # for k in range(1,self.nrows):
        #     for i in N:
        #         for j in N:
        #             if gx[i,j].x==k:
        #                 print(f"{k}: {i}->{j}")
        # print(self.order_rows, self.order_cols)

        transformed_matrix = self.matrix[self.order_rows,:][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix
    # def _solve_spp(self):

    #     # Modelcle
    #     model = gp.Model("tsp")
    #     model.setParam("Outputflag", self.output)
    #     model.setParam("TimeLimit", self.time_limit)
    #     #model.setParam("MIPFocus", 2)


    #     N=range(self.nrows)
    #     M=range(self.ncols)
    #     # Compute costs
    #     C1={}
    #     C2={ }
    #     for i in N:
    #         for j in N:
    #             C1[i, j] =  sum(abs(self.matrix[i, k] - self.matrix[j, k]) for k in M)
    #     for i in M:
    #         for j in M:
    #             C2[i, j] =  sum(abs(self.matrix[k, i] - self.matrix[k, j]) for k in N)


    #     # Variables
    #     x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
    #     y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        
    #     ux = model.addVars(self.nrows, lb=1, ub=self.nrows, name="ux")
    #     uy = model.addVars(self.ncols,lb=1, ub=self.ncols, name="sy")

    #     sx = model.addVars(self.nrows, vtype=GRB.BINARY,name="sx")
    #     sy = model.addVars(self.ncols, vtype=GRB.BINARY,name="sy")
    #     tx = model.addVars(self.nrows, vtype=GRB.BINARY,name="tx")
    #     ty = model.addVars(self.ncols, vtype=GRB.BINARY,name="ty")


    #     if self.coordinated_ordering and self.ncols==self.nrows:
    #         for i in N:
    #             for k in N:
    #                 model.addConstr(y[i,k]==x[i,k])
    #                 model.addConstr(gy[i,k]==gx[i,k])
    #             model.addConstr(zy[i]==zx[i])


    #     # Objective function
    #     model.setObjective(gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N) + gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M), GRB.MINIMIZE)


    #     # Constraints
    #     model.addConstr(sx[0] == 0)
    #     model.addConstr(gp.quicksum(sx[i] for i in N) == 1)
    #     model.addConstr(gp.quicksum(tx[i] for i in N) == 1)
    #     for i in N:
    #         model.addConstr(sx[i]+tx[i] <=1)
    #         model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) == 1-tx[i])
    #         x[i, i].ub=0
    #         model.addConstr(gp.quicksum(x[j,i] for j in N if i != j) == 1-sx[i])
    #         for j in N:
    #             model.addConstr(ux[i]-ux[j] + self.nrows*x[i,j] <= self.nrows-1)
    #             if j<i:
    #                 model.addConstr(x[i, j] + x[j, i] <= 1)
        
            
    #     if not self.coordinated_ordering:
    #         model.addConstr(sy[0] == 0)
    #         model.addConstr(gp.quicksum(sy[i] for i in M) == 1)
    #         model.addConstr(gp.quicksum(ty[i] for i in M) == 1)
    #         for i in M:
    #             model.addConstr(sy[i]+ty[i] <=1)
    #             model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) == 1-ty[i])
    #             model.addConstr(gp.quicksum(y[j, i] for j in M if i != j)  == 1-sy[i])
    #             model.addConstr(y[i, i] == 0)
    #             for j in M:
    #                 model.addConstr(uy[i]-uy[j] + self.ncols*y[i,j] <= self.ncols-1)
    #                 if j<i:
    #                     model.addConstr(y[i, j] + y[j, i] <= 1)


    #     # Solve model
    #     model._lproot=0
    #     model.optimize(self.lproot)
    #     #model.write("spp.sol")

    #     # Extract solution
    #     orderx=[]
    #     for i in N:
    #         if sx[i].x>0.5:
    #             orderx.append(i)
    #             i0=i
    #     while tx[i0].x<0.5:
    #         for j in N:
    #             if x[i0,j].x>0.5:
    #                 orderx.append(j)
    #                 i0=j
    #                 break

    #     if self.coordinated_ordering:
    #         ordery=orderx
    #     else:
    #         ordery=[]
    #         for i in M:
    #             if sy[i].x>0.5:
    #                 ordery.append(i)
    #                 i0=i
    #         while ty[i0].x<0.5:
    #             for j in M:
    #                 if y[i0,j].x>0.5:
    #                     ordery.append(j)
    #                     i0=j
    #                     break




    #     self.order_cols = ordery
    #     self.order_rows = orderx
        
    #     # for k in range(1,self.nrows):
    #     #     for i in N:
    #     #         for j in N:
    #     #             if gx[i,j].x==k:
    #     #                 print(f"{k}: {i}->{j}")
    #     # print(self.order_rows, self.order_cols)

    #     transformed_matrix = self.matrix[self.order_rows][:, self.order_cols]
    #     self.transformed_matrix=transformed_matrix
    #     self._extract_solution_info(model)

    #     return transformed_matrix
    
    # def _solve_spp8(self):

    #     # Modelcle
    #     model = gp.Model("tsp")
    #     model.setParam("Outputflag", self.output)
    #     model.setParam("TimeLimit", self.time_limit)
    #     #model.setParam("MIPFocus", 2)


    #     N=range(self.nrows)
    #     M=range(self.ncols)
    #     # Compute costs
    #     C1={}
    #     C2={ }
    #     for i in N:
    #         for j in N:
    #             C1[i, j] =  sum(abs(self.matrix[i, k] - self.matrix[j, k]) for k in M)
    #     for i in M:
    #         for j in M:
    #             C2[i, j] =  sum(abs(self.matrix[k, i] - self.matrix[k, j]) for k in N)


    #     # Variables
    #     x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
    #     y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        
    #     ux = model.addVars(self.nrows, lb=1, ub=self.nrows, name="ux")
    #     uy = model.addVars(self.ncols,lb=1, ub=self.ncols, name="sy")

    #     sx = model.addVars(self.nrows, vtype=GRB.BINARY,name="sx")
    #     sy = model.addVars(self.ncols, vtype=GRB.BINARY,name="sy")
    #     tx = model.addVars(self.nrows, vtype=GRB.BINARY,name="tx")
    #     ty = model.addVars(self.ncols, vtype=GRB.BINARY,name="ty")


    #     if self.coordinated_ordering:
    #         y=x
    #         sy=sx
    #         ty=tx
    #         uy=ux


    #     # Objective function
    #     C={}
    #     for i1 in N:
    #         for i2 in N:
    #             for j1 in M:
    #                 for j2 in M:
    #                     C[i1,i2,j1,j2]= np.abs(self.matrix[i1,j1]-self.matrix[i2,j2])+\
    #                     np.abs(self.matrix[i1,j2]-self.matrix[i2,j1])
                        

                        
    #    # Objective function
    #     obj1=gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N)
    #     obj2=gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M)
    #     obj12=gp.quicksum(C[i1,i2,j1,j2]*x[i1,i2]*y[j1,j2] for i1 in N for i2 in N for j1 in M for j2 in M)
    #     model.setObjective(obj1+obj2+obj12, GRB.MINIMIZE)


    #     # Constraints
    #     model.addConstr(sx[0] == 0)
    #     model.addConstr(gp.quicksum(sx[i] for i in N) == 1)
    #     model.addConstr(gp.quicksum(tx[i] for i in N) == 1)
    #     for i in N:
    #         model.addConstr(sx[i]+tx[i] <=1)
    #         model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) == 1-tx[i])
    #         x[i, i].ub=0
    #         model.addConstr(gp.quicksum(x[j,i] for j in N if i != j) == 1-sx[i])
    #         for j in N:
    #             model.addConstr(ux[i]-ux[j] + self.nrows*x[i,j] <= self.nrows-1)
    #             if j<i:
    #                 model.addConstr(x[i, j] + x[j, i] <= 1)
        
            
    #     if not self.coordinated_ordering:
    #         model.addConstr(sy[0] == 0)
    #         model.addConstr(gp.quicksum(sy[i] for i in M) == 1)
    #         model.addConstr(gp.quicksum(ty[i] for i in M) == 1)
    #         for i in M:
    #             model.addConstr(sy[i]+ty[i] <=1)
    #             model.addConstr(gp.quicksum(y[i, j] for j in M if j != i) == 1-ty[i])
    #             model.addConstr(gp.quicksum(y[j, i] for j in M if i != j)  == 1-sy[i])
    #             model.addConstr(y[i, i] == 0)
    #             for j in M:
    #                 model.addConstr(uy[i]-uy[j] + self.ncols*y[i,j] <= self.ncols-1)
    #                 if j<i:
    #                     model.addConstr(y[i, j] + y[j, i] <= 1)


    #     # Solve model
    #     model._lproot=0
    #     model.optimize(self.lproot)
    #     #model.write("spp.sol")

    #     # Extract solution
    #     orderx=[]
    #     for i in N:
    #         if sx[i].x>0.5:
    #             orderx.append(i)
    #             i0=i
    #     while tx[i0].x<0.5:
    #         for j in N:
    #             if x[i0,j].x>0.5:
    #                 orderx.append(j)
    #                 i0=j
    #                 break

    #     if self.coordinated_ordering:
    #         ordery=orderx
    #     else:
    #         ordery=[]
    #         for i in M:
    #             if sy[i].x>0.5:
    #                 ordery.append(i)
    #                 i0=i
    #         while ty[i0].x<0.5:
    #             for j in M:
    #                 if y[i0,j].x>0.5:
    #                     ordery.append(j)
    #                     i0=j
    #                     break




    #     self.order_cols = ordery
    #     self.order_rows = orderx
        
    #     # for k in range(1,self.nrows):
    #     #     for i in N:
    #     #         for j in N:
    #     #             if gx[i,j].x==k:
    #     #                 print(f"{k}: {i}->{j}")
    #     # print(self.order_rows, self.order_cols)

    #     transformed_matrix = self.matrix[self.order_rows][:, self.order_cols]
    #     self.transformed_matrix=transformed_matrix
    #     self._extract_solution_info(model)

    #     return transformed_matrix
    
    def _compute_error(self, A):
        
        error=np.zeros((self.nrows, self.ncols))
        for i in range(self.nrows):
            for j in range(self.ncols):
                error[i,j]=sum(np.abs(A[i,j]-A[a[0],a[1]]) for a in self.get_neighbors(i,j))

            
        return sum(sum(error)) if self.method in ["general", "s"] else 0.5*sum(sum(error))



    def _extract_solution_info(self, model):
        """Extracts relevant solution information after solving the model."""
        self.solution_info = {
            "file": self.file,
            "Formulation": self.method,
            "nrows": self.nrows,
            "ncols": self.ncols,
            "Sym": self.coordinated_ordering,
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
            # "Homogeneity Original": self.Homogeneity(self.matrix),
            # "Homogeneity Transformed": self.Homogeneity(self.transformed_matrix),
            "Order Rows": self.order_rows,
            "Order Cols": self.order_cols,
            # "Error Original": self._compute_error(self.matrix),
            # "Error Transformed": self._compute_error(self.transformed_matrix)
            #"TMatrix":  self.transformed_matrix
        }

    def get_neighbors(self, i, j):
        """Computes neighbors for a given (i, j) in the matrix."""

        if self.method=="VonNeumann":
            eps=1
        elif self.method=="Moore":
            eps=1.5
        else:
            eps=1

        neighbors=[]
        for i1 in range(self.nrows):
            for j1 in range(self.ncols):
                if 0.1 <= np.linalg.norm([i-i1, j-j1])<=eps:
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


    def plot_matrices(self, names=False):
        """Plots the original and transformed matrices side by side."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(18, 11))
        
        
        cmap = sns.color_palette("coolwarm", as_cmap=True)
        if names:
            axisAx=names
            axisAy=names
        else:
            axisAx=[i+1 for i in range(self.nrows)]
            axisAy=[i+1 for i in range(self.ncols)]
        vmin=self.matrix.min()
        vmax=self.matrix.max()
        heat1=sns.heatmap(self.matrix.T, annot=False, fmt=".2f", cmap=cmap, xticklabels=axisAx, yticklabels=axisAy, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[0], annot_kws={"size": 4}, cbar_kws={"shrink": 0.4})
        heat1.invert_yaxis()
        axes[0].set_xticklabels(axes[0].get_xticklabels(), fontsize=8)
        axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=8)
        axes[0].set_title(f"Original Matrix")
        axes[0].set_aspect("equal", adjustable="box")
        #axes[0].axis('equal')
        if names:
            axisBx=[names[i] for i in self.order_rows]
            axisBy=[names[i] for i in self.order_cols]
        else:
            axisBx=self.order_rows+np.ones(self.nrows, dtype=np.int8)
            axisBy=self.order_cols+np.ones(self.ncols, dtype=np.int8)

        heat2=sns.heatmap(self.transformed_matrix.T, annot=False, fmt=".2f", xticklabels=axisBx, yticklabels=axisBy, cmap=cmap, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[1], annot_kws={"size": 4}, cbar_kws={"shrink": 0.4})
        heat2.invert_yaxis()
        axes[1].set_title(f"Transformed Matrix")
        axes[1].set_xticklabels(axes[1].get_xticklabels(), fontsize=8)
        axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=8)
        axes[1].set_aspect("equal", adjustable="box")
        #axes[1].axis('equal')
        plt.show()
        #plt.savefig(f"{self.file[:-4]}_{self.method}_sym{self.coordinated_ordering}_eps{self.eps_neigh}.png", dpi=300, bbox_inches='tight')
        #plt.close()
    
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
    A=np.abs(p[:, np.newaxis] - p[np.newaxis, :])
    A=A-A.min()
    A=A/A.max()
    return np.round(A, decimals=2)

def generate_sqr_instances(n: int):
    p = np.random.uniform(0, 100, (n,2))
    A=np.sqrt(((p[:, np.newaxis, :] - p[np.newaxis, :, :]) ** 2).sum(axis=2))
    A=A-A.min()
    A=A/A.max()
    return np.round(A, decimals=2)

def generate_nonsqr_instances(n: int, m:int):
    p = np.random.uniform(0, 100, (n,2))
    q = np.random.uniform(0, 100, (m,2))
    A=np.sqrt(((p[:, np.newaxis, :] - q[np.newaxis, :, :]) ** 2).sum(axis=2))
    A=A-A.min()
    A=A/A.max()

    return np.round(A, decimals=2)


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

def compute_balanced_block_size(n, m, K):
    """
    Computes balanced block height and width for a given matrix size (n, m) and K value,
    ensuring the rectangles are as regular as possible and do not span the entire number of rows or columns.

    Parameters:
        n (int): Number of rows in the matrix.
        m (int): Number of columns in the matrix.
        K (int): Number of unique values (blocks).

    Returns:
        (int, int): Approximate block height and width.
    """
    total_blocks = K  # Aim for K different block values

    # Compute approximate block size
    block_h = max(1, n // int(np.sqrt(K)))  # Use square root scaling for balanced partitioning
    block_w = max(1, m // int(np.sqrt(K)))  # Ensure we get reasonable partitions

    # Ensure blocks do not cover the entire matrix
    block_h = min(block_h, n // 3)  # At most half of the rows
    block_w = min(block_w, m // 3)  # At most half of the columns

    return block_h, block_w

def generate_balanced_block_matrix(n, m, K):
    """
    Generates an n x m matrix with integer entries in [0, K], 
    initially structured as a block matrix with automatically determined block sizes,
    ensuring blocks do not span the entire number of rows or columns,
    then shuffles the rows and columns to disorder the structure.

    Parameters:
        n (int): Number of rows
        m (int): Number of columns
        K (int): Maximum integer value for entries (range: 0 to K)

    Returns:
        np.ndarray: Shuffled matrix
    """
    # Compute optimal block sizes
    block_h, block_w = compute_balanced_block_size(n, m, K)

    # Initialize an empty matrix
    matrix = np.zeros((n, m), dtype=int)

    # Assign values in computed blocks
    for i in range(0, n, block_h):
        for j in range(0, m, block_w):
            value = np.random.randint(0, K + 1)
            matrix[i:i + block_h, j:j + block_w] = value

    original=matrix
    # Shuffle rows and columns to disorder the matrix
    matrix = matrix[np.random.permutation(n), :]
    matrix = matrix[:, np.random.permutation(m)]

    return original, matrix

    



    