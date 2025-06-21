import numpy as np
import gurobipy as gp
from gurobipy import GRB
import seaborn as sns
import pandas as pd



def lproot(model, where):
        """
        Gurobi callback function to add cutting planes at MIPNODE and MIPSOL nodes.
        """
        if where == GRB.Callback.MIPNODE:  # Node solution
            nodecnt = model.cbGet(GRB.Callback.MIPNODE_NODCNT)
            if nodecnt==0:
                model._lp=model.cbGet(GRB.Callback.MIPNODE_OBJBND)

class MatrixSeriation:
    def __init__(self, matrix, file="test", method='general', time_limit=3600, coordinated_ordering=True, output=0, eps_neigh="Neumann"):
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
        self.coordinated_ordering = coordinated_ordering
        self.nrows, self.ncols = self.matrix.shape  # Rows and columns
        if self.nrows == self.ncols:
            self.square=True
        self.solution_info = {}
        self.output=output
        self.eps_neigh=eps_neigh
            
    def solve(self):
        """
        Solve the seriation problem using the selected method.
        """
        if self.method == 'PAM1':
            return self._solve_PAM1()
        elif self.method == 'PAM2':
            return self._solve_PAM2()
        elif self.method == 'HPM':
            if self.eps_neigh=="Neumann":
                return self._solve_HPM_VonNeumann()
            elif self.eps_neigh=="Moore":
                return self._solve_HPM_Moore()
            elif self.eps_neigh == "Cross":
                return self._solve_HPM_Cross()
            elif self.eps_neigh == "ME":
                return self._solve_HPM_ME()
            else:
                return "Non supported neighborhood"
        else:
            raise ValueError(f"Unknown method: {self.method}")
    
    
    
    def _solve_PAM1(self):

        """
        Solve the seriation problem using a linear ordering formulation.
        """
        model = gp.Model("Model_PAM1")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        

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
        model.optimize(lproot)

        if model.SolCount == 0:
            transformed_matrix = self.matrix.copy()
            self.order_rows = [i+1 for i in N]
            self.order_cols = [i+1 for i in M]
        else:
            transformed_matrix = np.empty((self.nrows,self.ncols))
            for k in N:
                for l in M:
                    transformed_matrix[k,l]=sum(self.matrix[i,j]*x[i,k].x*y[j,l].x for i in N for j in M)
            self.order_rows=[int(sum(i*x[i,k].x for i in N)) for k in N]
            self.order_cols=[int(sum(i*y[i,k].x for i in M)) for k in M]
    
        self.transformed_matrix=transformed_matrix
    

        self._extract_solution_info(model)

        return transformed_matrix
    
    
    def _solve_PAM2(self):

        # Model
        model = gp.Model("Model_PAM2")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)

        N=range(self.nrows)
        M=range(self.ncols)
        # Decision variables
        x = model.addVars(self.nrows, self.nrows, vtype=GRB.BINARY, name="x")
        y = model.addVars(self.ncols, self.ncols, vtype=GRB.BINARY, name="y")
        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])
        UB = self.matrix.max()
        s = model.addVars(self.nrows, self.ncols, ub=UB, name="s")
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
            s[k, l] >= UB*x[i, k] + gp.quicksum(self.matrix[i,j]*y[j, l] for j in M if self.matrix[i,j]>0.001) - UB
            for i in N for k in N for l in M
        )
        model.addConstrs(
            s[k, l] >= UB*y[j,l] + gp.quicksum(self.matrix[i,j]*x[i,k] for i in N if self.matrix[i,j]>0.001) - UB
            for j in M for k in N for l in M
        )
        model.addConstrs(
            s[k, l] <= UB*(1 - x[i,k]) + gp.quicksum(self.matrix[i,j]*y[j,l] for j in M if self.matrix[i,j]>0.001)
            for i in N for k in N for l in M
        )
        model.addConstrs(
            s[k, l] <= UB*(1 - y[j,l]) + gp.quicksum(self.matrix[i,j]*x[i,k] for i in N if self.matrix[i,j]>0.001)
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
        model.optimize(lproot)
        # Extract solution information
        #self._extract_solution_info(model)
        if model.Status==GRB.OPTIMAL or model.Status==GRB.SUBOPTIMAL:
        # transformed matrix
            transformed_matrix = np.empty((self.nrows,self.ncols))
            for k in N:
                for l in M:
                    transformed_matrix[k,l]=sum(self.matrix[i,j]*x[i,k].x*y[j,l].x for i in N for j in M)
        
            self.transformed_matrix=transformed_matrix
            self.order_rows=[int(sum((i)*x[i,k].x for i in N)) for k in N]
            self.order_cols=[int(sum((i)*y[i,k].x for i in M)) for k in M]
        else:
            transformed_matrix=self.matrix
            self.transformed_matrix=transformed_matrix
            self.order_rows=list(N)
            self.order_cols=list(M)


        self._extract_solution_info(model)

        return transformed_matrix



    def _solve_HPM_VonNeumann(self):

        # Modelcle
        model = gp.Model("HPM_VonNeumann")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("Outputflag", 1)
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
        model.optimize(lproot)

        if model.Status==GRB.OPTIMAL or model.Status==GRB.SUBOPTIMAL:
            #print("Obj: ", model.ObjVal)
            # Extract solution
            self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M) if zy[i].x < 0.5 else self.ncols-1 for i in M])
            self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N) if zx[i].x < 0.5 else self.nrows-1 for i in N])
            
            transformed_matrix = self.matrix[self.order_rows,:][:, self.order_cols]
            self.transformed_matrix=transformed_matrix
        else:
            transformed_matrix=self.matrix
            self.transformed_matrix=transformed_matrix
            self.order_rows=list(N)
            self.order_cols=list(M)


        self._extract_solution_info(model)

        return transformed_matrix
    

    def _solve_HPM_Moore(self):

        # Modelcle
        model = gp.Model("HPM_Moore")
        model.setParam("Outputflag", self.output)
        model.setParam("TimeLimit", self.time_limit)
        #model.setParam("MIPFocus", 3)




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

        

        gx = model.addVars(self.nrows, self.nrows, vtype=GRB.INTEGER, name="gx")
        gy = model.addVars(self.ncols, self.ncols, vtype=GRB.INTEGER, name="gy")

        zx = model.addVars(self.nrows, vtype=GRB.BINARY, name="zx")
        zy = model.addVars(self.ncols, vtype=GRB.BINARY, name="zy")

        if self.coordinated_ordering and self.ncols==self.nrows:
            for i in N:
                for k in N:
                    model.addConstr(y[i,k]==x[i,k])
                    model.addConstr(gy[i,k]==gx[i,k])
                model.addConstr(zy[i]==zx[i])



        theta = model.addVars(self.nrows, self.nrows, self.ncols, self.ncols, vtype=GRB.BINARY, name="theta")

        for i1 in N:
            for i2 in N:
                for j1 in M:
                    for j2 in M:
                        model.addConstr(theta[i1,i2,j1,j2] <=x[i1,i2])
                        model.addConstr(theta[i1,i2,j1,j2] <=y[j1,j2])
                        model.addConstr(theta[i1,i2,j1,j2] >= x[i1,i2] + y[j1,j2]-1)


        C={}
        for i1 in N:
            for i2 in N:
                for j1 in M:
                    for j2 in M:
                        C[i1,i2,j1,j2]= np.abs(self.matrix[i1,j1]-self.matrix[i2,j2])+ np.abs(self.matrix[i1,j2]-self.matrix[i2,j1])
                        

                        
       # Objective function
        obj1=gp.quicksum(C1[i, j]*x[i, j] for i in N for j in N)
        obj2=gp.quicksum(C2[i, j]*y[i, j] for i in M for j in M)
        obj12=gp.quicksum(C[i1,i2,j1,j2]*theta[i1,i2,j1,j2] for i1 in N for i2 in N for j1 in M for j2 in M)
        model.setObjective(obj1+obj2+obj12, GRB.MINIMIZE)


        # Constraints
        model.addConstr(zx[0] == 0)
        model.addConstr(gp.quicksum(zx[i] for i in N) == 1)
        for i in N:
            model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) + zx[i] == 1)
            model.addConstr(gp.quicksum(x[j,i] for j in N if i != j) <= 1)
            model.addConstr(gp.quicksum(gx[i, j] for j in N) >= gp.quicksum(gx[j, i] for j in N) - (self.nrows) * zx[i] + 1)
            for j in N:
                if j<i:
                    model.addConstr(x[i, j] + x[j, i] <= 1)
                model.addConstr(gx[i, j] <= (self.nrows-1) * x[i, j])
        
            
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
        model.optimize(lproot)

        if model.Status==GRB.OPTIMAL or model.Status==GRB.SUBOPTIMAL:
            # Extract solution
            self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M) if zy[i].x < 0.5 else self.ncols-1 for i in M])
            self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N) if zx[i].x < 0.5 else self.nrows-1 for i in N])
            
            transformed_matrix = self.matrix[self.order_rows,:][:, self.order_cols]
            self.transformed_matrix=transformed_matrix
        else:
            transformed_matrix=self.matrix
            self.transformed_matrix=transformed_matrix
            self.order_rows=list(N)
            self.order_cols=list(M)
        
        self._extract_solution_info(model)

        return transformed_matrix
    
    def _solve_HPM_Cross(self):

        # Modelcle
        model = gp.Model("HPM_cross")
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
        model.optimize(lproot)

        # Extract solution
        self.order_cols = np.argsort([sum(round(gy[i, j].x) for j in M) if zy[i].x < 0.5 else self.ncols for i in M])
        self.order_rows = np.argsort([sum(round(gx[i, j].x) for j in N) if zx[i].x < 0.5 else self.nrows for i in N])
        


        transformed_matrix = self.matrix[self.order_rows,:][:, self.order_cols]
        self.transformed_matrix=transformed_matrix
        self._extract_solution_info(model)

        return transformed_matrix
    
    def _solve_HPM_ME(self):

        # Modelcle
        model = gp.Model("HPM_ME")
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
        model.optimize(lproot)

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
            #"Memory": model.getAttr(GRB.Attr.MemUsed),
            "NumBinVars": model.NumBinVars,
            "NumCtrs": model.NumConstrs,
            "Order Rows": self.order_rows,
            "Order Cols": self.order_cols,
            "Homogeneity_NonTransf": self.Homogeneity(self.matrix),
            "Homogeneity_Transf": self.Homogeneity(self.transformed_matrix),
            "VNMeasure_NonTransf":  MeasuresExternal(self.matrix, range(self.nrows), range(self.ncols), index="Neumann"),#MeasureNT(index="Neumann"),
            "VNMeasure_Transf": MeasuresExternal(self.transformed_matrix, self.order_rows, self.order_cols, index="Neumann"),#self.MeasureT(index="Neumann"),
            "MooreMeasure_NonTransf": MeasuresExternal(self.matrix, range(self.nrows), range(self.ncols), index="Moore"),
            "MooreMeasure_Transf": MeasuresExternal(self.transformed_matrix, self.order_rows, self.order_cols, index="Moore"),
            "MEMeasure_NonTransf": MeasuresExternal(self.matrix, range(self.nrows), range(self.ncols), index="ME"),
            "MEMeasure_Transf": MeasuresExternal(self.transformed_matrix, self.order_rows, self.order_cols, index="ME"),
        }

    def get_neighbors(self, i, j):
        """Computes neighbors for a given (i, j) in the matrix."""

        if self.method=="Neumann":
            eps=1
        elif self.method=="Moore":
            eps=1.5
        else:
            eps=2

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
                HOM[i,j]=(1/len(self.get_neighbors(i,j)))*sum(abs(matrix[i,j]-matrix[i1,j1]) for (i1,j1) in self.get_neighbors(i,j))

        return (1/(self.nrows*self.ncols))*sum(sum(HOM))


    def plot_matrices(self, names=None, save=False):
        """Plots the original and transformed matrices side by side."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(20, 14)) 
        
        
        cmap = sns.color_palette("binary", as_cmap=True) #coolwarm Greys
        if names is not None and len(names) > 0:
            axisAx = names
            axisAy = names
        else:
            axisAx = [i+1 for i in range(self.nrows)]
            axisAy = [i+1 for i in range(self.ncols)]
        vmin=self.matrix.min()
        vmax=self.matrix.max()
        heat1=sns.heatmap(self.matrix.T, annot=False, fmt=".2f", cmap=cmap, xticklabels=axisAx, yticklabels=axisAy, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[0], annot_kws={"size": 4}, cbar_kws={"shrink": 0.4})
        heat1.invert_yaxis()
        axes[0].set_xticklabels(axes[0].get_xticklabels(), fontsize=8)
        axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=8)
        axes[0].set_aspect("equal", adjustable="box")
        extent = axes[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

        if save:
            fig.savefig(f"matrix1.png", bbox_inches=extent.expanded(1.1, 1.1), dpi=300)
        
        axes[0].set_title(f"Original Matrix")
        
        axisBx=self.order_rows+np.ones(self.nrows, dtype=np.int8)
        axisBy=self.order_cols+np.ones(self.ncols, dtype=np.int8)
        if names is not None and len(names) > 0:
            axisBx=[names[i] for i in self.order_rows]
            axisBy=[names[i] for i in self.order_cols]
            

        heat2=sns.heatmap(self.transformed_matrix.T, annot=False, fmt=".2f", xticklabels=axisBx, yticklabels=axisBy, cmap=cmap, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[1], annot_kws={"size": 4}, cbar_kws={"shrink": 0.4})
        heat2.invert_yaxis()
        
        axes[1].set_xticklabels(axes[1].get_xticklabels(), fontsize=8)
        axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=8)
        axes[1].set_aspect("equal", adjustable="box")
        extent = axes[1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        if save:
            fig.savefig(f"matrix2.png", bbox_inches=extent.expanded(1.1, 1.1), dpi=300)

        
        axes[1].set_title(f"Transformed Matrix")
        if save:
            fig.savefig(f"matrices.png", bbox_inches="tight", dpi=300)

        plt.show()
        #plt.savefig(f"{self.file[:-4]}_{self.method}_sym{self.coordinated_ordering}_eps{self.eps_neigh}.png", dpi=300, bbox_inches='tight')
        #plt.close()

    def MeasureT(self, index="Neumann"):

        #n, m = self.matrix.shape
        tmatrix=self.transformed_matrix
        return MeasuresExternal(tmatrix, self.order_rows, self.order_cols, index=index)
        
    
    def MeasureNT(self, index="Neumann"):

        #n, m = self.matrix.shape
        return MeasuresExternal(self.matrix, self.order_rows, self.order_cols, index=index)
        


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


def MeasuresExternal(matrix, order_rows, order_cols, index="Neumann"):

    n, m = matrix.shape
    #print("**", n, m, len(order_rows), len(order_cols))
    # print(f"n={n}, lenrows={len(order_rows)} m={m}, lencols={len(order_cols)}")
    # print(order_rows, order_cols)
    if n==len(order_rows) and m==len(order_cols):
        tmatrix=matrix[order_rows,:][:,order_cols]
        if index=="Neumann":
            eps=1
            neigh={}
            for i in range(n):
                for j in range(m):
                    neighbors=[]
                    for i1 in range(n):
                        for j1 in range(m):
                            if 0.1 <= np.linalg.norm(np.array([i-i1, j-j1]))<=eps:
                                neighbors.append((i1,j1))
                    neigh[i,j]=neighbors
            return sum(np.abs(tmatrix[i,j]-tmatrix[a[0],a[1]]) for i in range(n) for j in range(m) for a in neigh[i,j])
        elif index=="Moore":
            eps=1.5
            neigh={}
            for i in range(n):
                for j in range(m):
                    neighbors=[]
                    for i1 in range(n):
                        for j1 in range(m):
                            if 0.1 <= np.linalg.norm(np.array([i-i1, j-j1]))<=eps:
                                neighbors.append((i1,j1))
                    neigh[i,j]=neighbors
            return sum(np.abs(tmatrix[i,j]-tmatrix[a[0],a[1]]) for i in range(n) for j in range(m) for a in neigh[i,j])
        elif index=="ME":
            x={}
            for i in range(n):
                for j in range(m):
                    x[i,j]=tmatrix[i,j]
            for i in range(n):
                x[i,m]=0
                x[i,-1]=0
            for j in range(m):
                x[n,j]=0
                x[-1,j]=0
            return 0.5*sum(x[i,j]*(x[i,j+1] + x[i,j-1] + x[i+1,j]+ x[i-1,j]) for i in range(n) for j in range(m))
    else:
        return 0

def Plot_matrices_External(matrix, order_rows, order_cols, names=None):
        """Plots the original and transformed matrices side by side."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(20, 14)) #(18,11)
        
        nrows, ncols=matrix.shape
        cmap = sns.color_palette("binary", as_cmap=True) #coolwarm Greys
        axisAx = [i+1 for i in range(nrows)]
        axisAy = [i+1 for i in range(ncols)]
        if names is not None and len(names) > 0:
            axisAx = names
            axisAy = names
            
        vmin=matrix.min()
        vmax=matrix.max()
        heat1=sns.heatmap(matrix.T, annot=False, fmt=".2f", cmap=cmap, xticklabels=axisAx, yticklabels=axisAy, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[0], annot_kws={"size": 4}, cbar_kws={"shrink": 0.4})
        heat1.invert_yaxis()
        axes[0].set_xticklabels(axes[0].get_xticklabels(), fontsize=8)
        axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=8)
        axes[0].set_title(f"Original Matrix")
        axes[0].set_aspect("equal", adjustable="box")
        #axes[0].axis('equal')
        axisBx=order_rows+np.ones(nrows, dtype=np.int8)
        axisBy=order_cols+np.ones(ncols, dtype=np.int8)
        if names is not None and len(names) > 0:
            axisBx=[names[i] for i in order_rows]
            axisBy=[names[i] for i in order_cols]



        transformed_matrix=matrix[order_rows, :][:, order_cols]

        heat2=sns.heatmap(transformed_matrix.T, annot=False, fmt=".2f", xticklabels=axisBx, yticklabels=axisBy, cmap=cmap, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes[1], annot_kws={"size": 4}, cbar_kws={"shrink": 0.4})
        heat2.invert_yaxis()
        axes[1].set_title(f"Transformed Matrix")
        axes[1].set_xticklabels(axes[1].get_xticklabels(), fontsize=8)
        axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=8)
        axes[1].set_aspect("equal", adjustable="box")
        plt.show()


def Plot_matrices_External_OnlyOne(matrix, order_rows, order_cols, names=None, namefile=""):
        """Plots the original and transformed matrices side by side."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 1, figsize=(20, 14)) #(18,11)
        
        nrows, ncols=matrix.shape
        cmap = sns.color_palette("binary", as_cmap=True) #coolwarm Greys
        vmin=matrix.min()
        vmax=matrix.max()
        axisBx=order_rows+np.ones(nrows, dtype=np.int8)
        axisBy=order_cols+np.ones(ncols, dtype=np.int8)
        if names is not None and len(names) > 0:
            axisBx=[names[i] for i in order_rows]
            axisBy=[names[i] for i in order_cols]
            
        transformed_matrix=matrix[order_rows, :][:, order_cols]

        heat2=sns.heatmap(transformed_matrix.T, annot=False, fmt=".2f", xticklabels=axisBx, yticklabels=axisBy, cmap=cmap, linewidths=0.5, vmin=vmin, vmax=vmax, cbar=True, center=0.5, ax=axes, annot_kws={"size": 4}, cbar_kws={"shrink": 0.4})
        heat2.invert_yaxis()
        #axes.set_title(f"Transformed Matrix")
        axes.set_xticklabels(axes.get_xticklabels(), fontsize=8)
        axes.set_yticklabels(axes.get_yticklabels(), fontsize=8)
        axes.set_aspect("equal", adjustable="box")
        #axes[1].axis('equal')
        
        if namefile:
            plt.savefig(namefile, dpi=300, bbox_inches='tight')
        plt.show()
        #plt.close()

def HomogeneityExternal(matrix, index="Neumann"):

    nrows, ncols=matrix.shape
    N=range(nrows)
    M=range(ncols)
    HOM=np.empty((nrows,ncols))


    if index=="Neumann":
            eps=1
            neigh={}
            for i in range(nrows):
                for j in range(ncols):
                    neighbors=[]
                    for i1 in range(nrows):
                        for j1 in range(ncols):
                            if 0.1 <= np.linalg.norm(np.array([i-i1, j-j1]))<=eps:
                                neighbors.append((i1,j1))
                    neigh[i,j]=neighbors
    elif index=="Moore":
        eps=1.5
        neigh={}
        for i in range(nrows):
            for j in range(ncols):
                neighbors=[]
                for i1 in range(nrows):
                    for j1 in range(ncols):
                        if 0.1 <= np.linalg.norm(np.array([i-i1, j-j1]))<=eps:
                            neighbors.append((i1,j1))
                neigh[i,j]=neighbors
    else:
        return 0

    for i in N:
        for j in M:
            HOM[i,j]=(1/len(neigh[i,j]))*sum(abs(matrix[i,j]-matrix[i1,j1]) for (i1,j1) in neigh[i,j])

    return (1/(nrows*ncols))*HOM.sum()

def random_partition(n, m):
    cuts = np.sort(np.random.randint(1, n, m - 1))  # Generate m-1 random cut points
    parts = np.diff([0] + cuts.tolist() + [n])  # Compute differences between cut points
    return parts

