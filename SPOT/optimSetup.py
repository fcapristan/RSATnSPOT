def IPOPTrun(vars0,mission,includeDrag):
    import ipopt
    import constraints as ct
    import costFunction as CF
    nvar = len(vars0)


    con1= ct.eval_g(vars0,mission)
    g_L,g_U = ct.bounds(vars0,mission)
    ncon = len(g_L)

    jac = ct.eval_jac_g(vars0,mission,'1d')

    #print np.shape(jac),ncon,len(vars0)
    check  = CF.fprimeObjectiveMass(vars0,mission)
    x_L,x_U = ct.getXLXU(vars0,mission)

    # start IPOPT classess
    class trajectory(object):
        def __init__(self,mission,includeDrag):
            self._mission = mission
            self._includeDrag = includeDrag
        def objective(self,x):
            ret = CF.ObjectiveMass(x,self._mission)
            return ret
        def gradient(self,x):
            ret = CF.fprimeObjectiveMass(x,self._mission)
            return ret
        def constraints(self,x):
            ret = ct.eval_g(x,self._mission,self._includeDrag)
            return ret
        def jacobian(self,x):
            ret = ct.eval_jac_g(x,self._mission,'1d',self._includeDrag)
            return ret
    nlp = ipopt.problem(
                        n=len(vars0),
                        m=ncon,
                        problem_obj=trajectory(mission,includeDrag),
                        lb=x_L,
                        ub=x_U,
                        cl=g_L,
                        cu=g_U
                        )

    #
    # Set solver options
    #
    nlp.addOption('mu_strategy', 'adaptive')
    nlp.addOption('tol', 1e-5)
    nlp.addOption('max_iter',1000)
    nlp.addOption('acceptable_constr_viol_tol',5e-6)
    nlp.addOption('constr_viol_tol',1e-5)
    #nlp.addOption('derivative_test','first-order')
    nlp.addOption('acceptable_tol',5e-6)
    nlp.addOption('expect_infeasible_problem','yes')
    '''
    #
    # Set solver options
    #
    nlp.addOption('mu_strategy', 'adaptive')
    nlp.addOption('tol', 5e-6)
    nlp.addOption('max_iter',50000)
    nlp.addOption('acceptable_constr_viol_tol',5e-7)
    nlp.addOption('constr_viol_tol',5e-7)
    nlp.addOption('derivative_test','first-order')
    nlp.addOption('acceptable_tol',5e-7)
    nlp.addOption('expect_infeasible_problem','yes')
    '''
    x,info = nlp.solve(vars0)
    return x,info


def SNOPTrun(vars0,mission,includeDrag,flagFORCE=1):
    import pyOpt
    import constraints
    import costFunction
    import numpy as np
            
    # THESE FUNCTIONS ARE USED BY SNOPT
    # Define the functions SNOPT optimizer will call


    def objectiveFunction(inVars,mission,includeDrag):
        x=costFunction.ObjectiveMass(inVars,mission)
        eq=constraints.equality(inVars,mission,'real',includeDrag)
        ineq=constraints.inequality(inVars,mission,'real')
        
        g = np.concatenate((eq,ineq),1)
        fail = 0
        return x,g,fail

    def sensitivityFunction(inVars,f,g,mission,includeDrag):
        x = costFunction.fprimeObjectiveMass(inVars,mission)
        eq = constraints.fprimeequality(inVars,mission,'2d',includeDrag)
        ineq = constraints.fprimeinequality(inVars,mission,'2d')
        
        g = np.concatenate((eq,ineq),0)
        fail = 0
        return x,g,fail


    numEquality = len(constraints.equality(vars0,mission))
    numInequality = len(constraints.inequality(vars0,mission))

    # Find the upper and lower bounds
    #boundsCase = constraints.bounds(mission)
    lb,ub = constraints.getXLXU(vars0,mission)


    #TJC
    opt_prob = pyOpt.Optimization('Trajectory Optimization',lambda x: objectiveFunction(x,mission,includeDrag))
    opt_prob.addObj('Objective Mass')

    # Specify all of the variables
    #print 'Setting up variables in a hackish way.  MUST CHANGE!!!'
    for curVar in range(len(vars0)):
     opt_prob.addVar('var'+str(curVar), 'c', value=vars0[curVar], lower=lb[curVar], upper=ub[curVar])

    # Now add in equality constraints
    for curCon in range(numEquality):
     opt_prob.addCon('g' + str(curCon), 'e')

    # Now add in inequality constraints
    for curCon in range(numEquality,numEquality + numInequality):
     opt_prob.addCon('g' + str(curCon), 'i')

    # Confirm that everything is correct
    #print opt_prob

    # Set up the optimizer
    snopt = pyOpt.pySNOPT.SNOPT()
    snopt.setOption('Major feasibility tolerance',value=5e-6)
    snopt.setOption('Major optimality tolerance',value=1e-5)
    snopt.setOption('Minor feasibility tolerance',value=5e-6)
    snopt.setOption('Major iterations limit',500)
    print 'Using SNOPT'

    # Optimize and save results
    sens2 = lambda x,f,g:sensitivityFunction(x,f,g,mission,includeDrag)

    # by default will try complex step first...if fails...then finite diference
    exitVal = snopt(opt_prob,sens_type= sens2)
    
    infoOpt = exitVal[2]['text']
    if infoOpt!='finished successfully' and flagFORCE==1:
       print 'Failed to finish successfully with CS .... trying FD'
       exitVal = snopt(opt_prob,sens_type= 'FD')

    return exitVal


def IPOPTrunINIT(vars0,mission,includeDrag):
    import ipopt
    import constraints as ct
    import costFunction as CF
    nvar = len(vars0)
    
    
    con1= ct.eval_g(vars0,mission)
    g_L,g_U = ct.bounds(vars0,mission)
    ncon = len(g_L)
    
    jac = ct.eval_jac_g(vars0,mission,'1d')
    
    #print np.shape(jac),ncon,len(vars0)
    check  = CF.fprimeObjectiveMassINIT(vars0,mission)
    x_L,x_U = ct.getXLXU(vars0,mission)
    
    # start IPOPT classess
    class trajectory(object):
        def __init__(self,mission,includeDrag):
            self._mission = mission
            self._includeDrag = includeDrag
        def objective(self,x):
            ret = CF.ObjectiveMassINIT(x,self._mission)
            return ret
        def gradient(self,x):
            ret = CF.fprimeObjectiveMassINIT(x,self._mission)
            return ret
        def constraints(self,x):
            ret = ct.eval_g(x,self._mission,self._includeDrag)
            return ret
        def jacobian(self,x):
            ret = ct.eval_jac_g(x,self._mission,'1d',self._includeDrag)
            return ret
    nlp = ipopt.problem(
                        n=len(vars0),
                        m=ncon,
                        problem_obj=trajectory(mission,includeDrag),
                        lb=x_L,
                        ub=x_U,
                        cl=g_L,
                        cu=g_U
                        )
    
    
    #
    # Set solver options
    #
    nlp.addOption('mu_strategy', 'adaptive')
    nlp.addOption('tol', 5e-6)
    nlp.addOption('max_iter',1000)
    nlp.addOption('acceptable_constr_viol_tol',5e-6)
    nlp.addOption('constr_viol_tol',5e-6)
    #nlp.addOption('derivative_test','first-order')
    nlp.addOption('acceptable_tol',5e-6)
    nlp.addOption('expect_infeasible_problem','yes')
    
    x,info = nlp.solve(vars0)
    return x,info







def SNOPTrunINIT(vars0,mission,includeDrag):
    import pyOpt
    import constraints
    import costFunction
    import numpy as np
    
    # THESE FUNCTIONS ARE USED BY SNOPT
    # Define the functions SNOPT optimizer will call
    
    
    def objectiveFunction(inVars,mission,includeDrag):
        x=costFunction.ObjectiveMassINIT(inVars,mission)
        eq=constraints.equality(inVars,mission,'real',includeDrag)
        ineq=constraints.inequality(inVars,mission,'real')
        
        g = np.concatenate((eq,ineq),1)
        fail = 0
        return x,g,fail
    
    def sensitivityFunction(inVars,f,g,mission,includeDrag):
        x = costFunction.fprimeObjectiveMassINIT(inVars,mission)
        eq = constraints.fprimeequality(inVars,mission,'2d',includeDrag)
        ineq = constraints.fprimeinequality(inVars,mission,'2d')
        
        g = np.concatenate((eq,ineq),0)
        fail = 0
        return x,g,fail
    
    
    numEquality = len(constraints.equality(vars0,mission))
    numInequality = len(constraints.inequality(vars0,mission))
    
    # Find the upper and lower bounds
    #boundsCase = constraints.bounds(mission)
    lb,ub = constraints.getXLXU(vars0,mission)
    
    
    #TJC
    opt_prob = pyOpt.Optimization('Trajectory Optimization',lambda x: objectiveFunction(x,mission,includeDrag))
    opt_prob.addObj('Objective Mass')
    
    # Specify all of the variables
    #print 'Setting up variables in a hackish way.  MUST CHANGE!!!'
    for curVar in range(len(vars0)):
        opt_prob.addVar('var'+str(curVar), 'c', value=vars0[curVar], lower=lb[curVar], upper=ub[curVar])
    
    # Now add in equality constraints
    for curCon in range(numEquality):
        opt_prob.addCon('g' + str(curCon), 'e')
    
    # Now add in inequality constraints
    for curCon in range(numEquality,numEquality + numInequality):
        opt_prob.addCon('g' + str(curCon), 'i')
    
    # Confirm that everything is correct
    #print opt_prob
    
    # Set up the optimizer
    snopt = pyOpt.pySNOPT.SNOPT()
    snopt.setOption('Major feasibility tolerance',value=5e-7)
    snopt.setOption('Major optimality tolerance',value=1e-5)
    snopt.setOption('Minor feasibility tolerance',value=5e-7)
    print 'Using SNOPT'

    # Optimize and save results
    sens2 = lambda x,f,g:sensitivityFunction(x,f,g,mission,includeDrag)
    exitVal = snopt(opt_prob,sens_type= sens2)
    #exitVal = snopt(opt_prob,sens_type= 'FD')
    

    return exitVal























