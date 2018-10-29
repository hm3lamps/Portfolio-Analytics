

import pandas as pd
import numpy as np

def readExcel(filename):
    xl = pd.ExcelFile(filename)
    df = xl.parse("Excess Returns")
    return df

def blackLitterman(MarketWeights, Q, tauOmega, casename, views):
    print("Processing--" + casename)
    df = readExcel("BL+-returndata.xlsx")
    df = df.iloc[:,2:5]
    covMat = np.cov(df, rowvar=False)
    RiskAversion = 3.0
    tau = 0.1
    
    MarketVar = float (np.matmul(np.matrix(MarketWeights).T,np.matmul ( covMat , MarketWeights)))
    MarketExpExRet = float ( RiskAversion * MarketVar )
    StdDev =  MarketVar **(1/2)
    sharpeRat = MarketExpExRet / StdDev
    priorReturns = RiskAversion * np.matmul( covMat , MarketWeights)
    covtau = tau * covMat
   
    omega = tauOmega*np.matmul(views,np.matmul(covMat , views.T) )
    diagOmega = np.diag(np.diag (omega))
    PriorPrecision = np.matmul( views, np.matmul( tau * covMat , views.T) )
    PosteriorReturns = priorReturns + np.matmul(np.matmul(covtau,
    np.matmul(views.T, np.linalg.inv(PriorPrecision+omega))),Q-np.matmul(
      views, priorReturns))
    PosteriorReturnDist = covMat + (covtau - np.matmul(np.matmul(np.matmul(
    covtau,views.T),np.linalg.inv(PriorPrecision+omega)),np.matmul(views,
    covtau)))   
    unconstrainedOpt = np.matmul( np.matrix(PosteriorReturns).T ,np.linalg.inv(RiskAversion * PosteriorReturnDist))
    sum = unconstrainedOpt.sum(dtype = 'float')
    OptWeight = np.matrix( unconstrainedOpt / sum ).T
    ExpectedReturn = float( np.matmul( np.matrix(PosteriorReturns).T, OptWeight) )
    variance =float(np.matmul( np.matmul (OptWeight.T , PosteriorReturnDist),OptWeight))
    StdDevOpt = variance**(1/2)
    SharpeOpt = ExpectedReturn/StdDevOpt
   
    rowheader = ['US Equity', 'Foreign Equity', 'Emerging Equity']
    covstring = "Covariance Matrix\n"
    covtaustring = "Covariance * tau\n"
    postretdiststring  = "Posterior Retrun Distributiom\n"
    covtaustring += ","
    covstring +="," 
    postretdiststring +=","

    for i in range (len(covMat)):
        covstring += "%s," % rowheader[i]    #creating col headers
        covtaustring += "%s," % rowheader[i]  
        postretdiststring += "%s," % rowheader[i]  
    covstring +="\n"
    covtaustring +="\n"
    postretdiststring +="\n"

    for i in range(len(covMat)):
        covstring += "%s," % rowheader[i] #creating row headers
        covtaustring += "%s," %rowheader[i]
        postretdiststring += "%s," %rowheader[i]

        for j in range(len(covMat[0])):
            covstring += "%g," % covMat[i][j]
            covtaustring += "%g," %covtau[i][j]
            postretdiststring += "%g," %PosteriorReturnDist[i][j]
        covstring += "\n"
        covtaustring+= "\n"
        postretdiststring+= "\n"
    covstring += "\n"
    covtaustring += "\n"
    postretdiststring+= "\n"
    covstring +="\n"
    covtaustring += "\n"
    postretdiststring+= "\n"

       # market weights, Prior Returns and Posterior Returns -> making them 
    Tmarketweights = np.array(MarketWeights).T
    TPriorRet   =  np.array(priorReturns).T
    TPostRet =  np.array(PosteriorReturns).T
    Toptweights = OptWeight.T
    marketstring = "Market Weights\n" 
    PriorRetstring = "Prior Returns\n"
    PostRetstring = "Posterior Returns\n"
    optweightstring = "Unconstrained Optimization using Matrices\n"
    constrainedstring = "Optimized Weight constrained to 1\n"

    for i in range (len(MarketWeights)):
        marketstring += "%s," % rowheader[i]    #creating col headers
        PriorRetstring += "%s," % rowheader[i] 
        PostRetstring += "%s," % rowheader[i] 
        optweightstring +="%s," % rowheader[i] 
        constrainedstring +="%s," % rowheader[i]

    marketstring +="\n"
    PriorRetstring +=  "\n "
    PostRetstring += "\n"
    optweightstring += "\n"
    constrainedstring += "\n"


    for i in range(len(Tmarketweights)):
        for j in range(len(Tmarketweights[0])):
            marketstring += "%g," % Tmarketweights[i][j]
            PriorRetstring += "%g," % TPriorRet[i][j]
            PostRetstring += "%g," % TPostRet[i][j]
            optweightstring += '%g,' % unconstrainedOpt[i, j]
            constrainedstring += '%g,' %Toptweights[i, j]
        marketstring += "\n"
        PriorRetstring +="\n"
        PostRetstring +="\n"
        optweightstring += "\n"
        constrainedstring += "\n"    

    marketstring += "\n"
    PriorRetstring+="\n"
    PostRetstring+="\n"
    optweightstring += "\n"
    constrainedstring += "\n"

    marketstring +="\n"
    PriorRetstring+="\n"
    PostRetstring+="\n"
    optweightstring += "\n"
    constrainedstring += "\n"

        #printing views 2*3 array 
    viewstring = "P,"
    viewrow =['view1', 'view2']
    for i in range(len(MarketWeights)):
        viewstring   += '%s,' %rowheader[i]
    viewstring += '\n'

    for i in range(len(views)):
        viewstring += '%s,' %viewrow[i]
        for j in range(len(MarketWeights)):
            viewstring += '%g,' %views[i][j]
        viewstring +='\n'
    viewstring +='\n'
    viewstring +='\n'

       #printing Q 1*2
    Qstring = ",Q\n"
    for i in range(len(Q)):
        Qstring += '%s,' %rowheader[i]
        for j in range(len(Q[0])):
            Qstring += '%g,' %Q[i][j]
        Qstring +='\n'
    Qstring +='\n'
    Qstring +='\n'

    omegastring = "Omega\n"
    diagomegastring = "Diagonal Omega\n"
    priorprecstring = "Prior Precision of Views\n"
    for i in range(len(omega)):
        omegastring += '%s,' %rowheader[i]
        diagomegastring += '%s,' %rowheader[i]
        priorprecstring += '%s,' %rowheader[i]
        for j in range(len(omega[0])):
            omegastring += '%g,' % omega[i][j]
            diagomegastring += '%g,' %diagOmega[i][j] 
            priorprecstring += '%g,' %PriorPrecision[i][j]
        omegastring +='\n'
        diagomegastring += '\n'
        priorprecstring += '\n'     
    omegastring +='\n'
    diagomegastring += '\n'
    priorprecstring += '\n'

    omegastring +='\n'
    diagomegastring += '\n'
    priorprecstring += '\n'
   
   
    filename = "BL" + casename
    with open(filename + ".csv", "a") as f:
        #f.write('Hardik Maheshwari\nGTID: 903216069\nBlack Litterman Assignment\n')
        #f.write('\n'*3)  
        #f.write('PRIOR\n')
        #f.write(covstring)

        f.write('Risk Aversion,')
        f.write(str(RiskAversion))
        f.write('\n')

        f.write('Tau,')
        f.write(str(tau))

        f.write('\nMarket Variance,')
        f.write(str(MarketVar))

        f.write('\nMarket Exp. Excess Return,')
        f.write(str(MarketExpExRet))

        f.write('\nStandard Deviation,')
        f.write(str(StdDev))

        f.write('\nSharpe Ratio,')
        f.write(str(sharpeRat))


        f.write('\n'*3)
        f.write(marketstring)
        f.write(PriorRetstring)
        f.write (covtaustring)

        f.write('\n'*2)
        f.write('VIEWS\n')
        f.write(viewstring)  
        f.write(Qstring)
        f.write(omegastring)
        f.write(diagomegastring)
        f.write(priorprecstring)

        f.write ('\n'*2)
        f.write('POSTERIOR ESTIMATES\n')
        f.write(PostRetstring)
        f.write(postretdiststring)

        f.write('\n'*2)
        f.write('Portfolio Optimization\n')
        f.write(optweightstring)
        f.write(constrainedstring)
        f.write('\nExpected return,')
        f.write(str(ExpectedReturn))   
        f.write('\nOpt Variance,')
        f.write(str(variance))
        f.write('\nOpt Standard Dev,')
        f.write(str(StdDevOpt))
        f.write('\nOpt Sharpe ratio,')
        f.write(str(SharpeOpt))

        f.close()
        print("Finshed Processing--" + casename)
    
if __name__ == '__main__':
    
    
    MarketWeights = [[0.5],[0.4],[0.1]]
    Q = [ [0.015], [0.030]]
    views = np.array([ [1, 0, 0], [0, 1, -1] ])
    
    #case1: 
    tauOmega = 0.1
    blackLitterman(MarketWeights, Q, tauOmega, "case1", views)
    
    #case2: 
    tauOmega = 0.01
    blackLitterman(MarketWeights, Q, tauOmega, "case2", views)
    
    #case3: 
    tauOmega = 0.1
    Q = [ [0.02], [0.015]]
    views = np.array([ [1, -1, 0], [0, 0, 1] ])
    blackLitterman(MarketWeights, Q, tauOmega, "case3", views)
    

    

