""" ANOVA GAGE R&R PYTHON MODULE
# Description:
    This is python(3.6) implementation of gageRRDesign(in R).
# Dependencies: Numpy, Pandas, Matplotlib
# Author: Shin-Fu (Kelvin) Wu
# Date: 2017/06/08
# Reference:
    * https://www.rdocumentation.org/packages/qualityTools/versions/1.31.1/topics/gageRRDesign
    * https://www.spcforexcel.com/knowledge/measurement-systems-analysis/anova-gage-rr-part-1
    * https://www.spcforexcel.com/knowledge/control-chart-basics/control-limits
"""
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class AnovaGRR():
    def __init__(self, measurement, sig=6, tolerance=None):
        """ ANOVA Gage R&R
        Parameters
        ----------
        measurement: dict{string appraiser: dict{int part : list trials}}
            The dictionary contains all the measurement results.
        sig: int/float, optional
            Magnification of standard deviation.
        tolerance: int/float, optional
            Tolerance of the measurement. It is for CL and P/T ratio.
        Attributes
        ----------
        k: int
            Number of appraisers (operators).
        r: int
            Number of trials (replications). 
        n: int
            Number of parts.
        ANOVAtable: pandas.DataFrame
            The table measures the variability from five sources (operator, part, interaction between operator and part, equipment, total)
        GRRtable: pandas.DataFrame
            Output results of the Gage R&R.
        """
        self.measurement = measurement
        self.sig = sig
        self.tolerance = tolerance
        self.__ANOVA()
        self.__GRR()
    
    def __ANOVA(self):
        self.k = len(self.measurement.keys())
        self.r = len(list(list(self.measurement.values())[0].values())[0])
        self.n = len(list(list(self.measurement.values())[0].values()))
        
        total = np.array([])
        for i in range(len(self.measurement.keys())):
            if i == 0:
                total = np.array(list(list(self.measurement.values())[i].values()))
            else:
                total = np.vstack((total, np.array(list(list(self.measurement.values())[i].values()))))
        Xbar = np.mean(total)
        SSt = 0
        for i in range(self.k):
            for j in range(self.n):
                for m in range(self.r):
                    app = self.measurement[list(self.measurement.keys())[i]]
                    Xijm = app[list(app.keys())[j]][m]
                    SSt += (Xijm - Xbar)**2
        SSo = 0
        for i in range(self.k):
            SSo += self.n*self.r*(np.mean(list(self.measurement[list(self.measurement.keys())[i]].values())) - Xbar)**2
        
        SSp = 0
        for j in range(self.n):
            X_j_ = []
            for i in range(self.k):
                app = self.measurement[list(self.measurement.keys())[i]]
                X_j_.append(app[list(app.keys())[j]])      
            SSp += self.k*self.r* (np.mean(X_j_) - Xbar)**2
        
        SSe = 0
        for i in range(self.k):
            for j in range(self.n):
                for m in range(self.r):
                    app = self.measurement[list(self.measurement.keys())[i]]
                    Xijbar = np.mean(app[list(app.keys())[j]])
                    SSe += (app[list(app.keys())[j]][m] - Xijbar)**2
        
        SSop = SSt - (SSo + SSp + SSe)
        
        DFo = self.k - 1
        DFp = self.n - 1
        DFe = self.n*self.k*(self.r-1)
        DFt = self.n*self.k*self.r - 1
        DFop = DFt - (DFo + DFp + DFe) 
            
        MSo = SSo/DFo
        MSp = SSp/DFp
        MSop = SSop/DFop
        MSe = SSe/DFe    
        table = np.array([
                [DFo, SSo, MSo, MSo/MSop],
                [DFp, SSp, MSp, MSp/MSop],
                [DFop, SSop, MSop, MSop/MSe],
                [DFe, SSe, SSe/DFe, None],
                [DFt, SSt, None, None]
                ])
        self.ANOVAtable = pd.DataFrame(table, columns=['DF', 'SS', 'MS', 'F'], index=['Operator','Part','Operator by Part','Equipment','Total'])

    def __GRR(self):
        MSo = self.ANOVAtable['MS']['Operator']
        MSp = self.ANOVAtable['MS']['Part']
        MSop = self.ANOVAtable['MS']['Operator by Part']
        MSe = self.ANOVAtable['MS']['Equipment']
        
        S2re = max(MSe, 0)
        S2op = max((MSop - S2re) / self.r, 0)
        S2p = max((MSp - MSop) / (self.k * self.r), 0)
        S2o = max((MSo - MSop) / (self.n * self.r), 0)
        
        GRR2 = S2re + S2o
        TV = S2re + S2op + S2o + S2p
        
        VarComp = GRR2
        VarCompContrib = GRR2/TV
        Stdev = np.sqrt(GRR2)
        StudyVar = self.sig*Stdev
        TotalVarStudyVar = self.sig*np.sqrt(TV)
        StudyVarContrib = StudyVar/TotalVarStudyVar
        
        if self.tolerance == 0 or self.tolerance == None:
            table = \
            [
             [VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib],
             [S2re, S2re/TV, np.sqrt(S2re), self.sig*np.sqrt(S2re), self.sig*np.sqrt(S2re)/TotalVarStudyVar],
             [S2o+S2op, (S2o+S2op)/TV, np.sqrt(S2o+S2op), self.sig*np.sqrt(S2o+S2op), self.sig*np.sqrt(S2o+S2op)/TotalVarStudyVar],
             [S2o, S2o/TV, np.sqrt(S2o), self.sig*np.sqrt(S2o), self.sig*np.sqrt(S2o)/TotalVarStudyVar],
             [S2op, S2op/TV, np.sqrt(S2op), self.sig*np.sqrt(S2op), self.sig*np.sqrt(S2op)/TotalVarStudyVar],
             [S2p, S2p/TV, np.sqrt(S2p), self.sig*np.sqrt(S2p), self.sig*np.sqrt(S2p)/TotalVarStudyVar],
             [TV, TV/TV, np.sqrt(TV), self.sig*np.sqrt(TV), self.sig*np.sqrt(TV)/TotalVarStudyVar]
            ]
            self.GRRtable = pd.DataFrame(table, columns=['VarComp','VarCompContrib','Stdev','StudyVar','StudyVarContrib'],\
                                         index=['TotalRR','Repeat.','Reprod.','    Operator', '    Op.:Part', 'Part to part', 'TotalVar'])
        else:
            PTratio = StudyVar/self.tolerance
            
            table = \
            [
             [VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, PTratio],
             [S2re, S2re/TV, np.sqrt(S2re), self.sig*np.sqrt(S2re), self.sig*np.sqrt(S2re)/TotalVarStudyVar, self.sig*np.sqrt(S2re)/self.tolerance],
             [S2o+S2op, (S2o+S2op)/TV, np.sqrt(S2o+S2op), self.sig*np.sqrt(S2o+S2op), self.sig*np.sqrt(S2o+S2op)/TotalVarStudyVar, self.sig*np.sqrt(S2o+S2op)/self.tolerance],
             [S2o, S2o/TV, np.sqrt(S2o), self.sig*np.sqrt(S2o), self.sig*np.sqrt(S2o)/TotalVarStudyVar, self.sig*np.sqrt(S2o)/self.tolerance],
             [S2op, S2op/TV, np.sqrt(S2op), self.sig*np.sqrt(S2op), self.sig*np.sqrt(S2op)/TotalVarStudyVar, self.sig*np.sqrt(S2op)/self.tolerance],
             [S2p, S2p/TV, np.sqrt(S2p), self.sig*np.sqrt(S2p), self.sig*np.sqrt(S2p)/TotalVarStudyVar, self.sig*np.sqrt(S2p)/self.tolerance],
             [TV, TV/TV, np.sqrt(TV), self.sig*np.sqrt(TV), self.sig*np.sqrt(TV)/TotalVarStudyVar, self.sig*np.sqrt(TV)/self.tolerance]
            ]
            self.GRRtable = pd.DataFrame(table, columns=['VarComp','VarCompContrib','Stdev','StudyVar','StudyVarContrib','P/T ratio'],\
                                         index=['TotalRR','Repeat.','Reprod.','    Operator', '    Op.:Part', 'Part to part', 'TotalVar'])
        
    def PrintResults(self):
        """ Print Results """
        print('===')
        print('ANOVA Table')
        print(self.ANOVAtable)
        print('===')
        print('Gage R&R')
        print(self.GRRtable)
        print('===')
    
    def PlotCharts(self, size=(10, 10), filename=None):
        """ Plot the Control Charts
        Parameters
        ----------
        size: tuple(width, height), optional
            Size of the figure.
        filename: string, optional
            the file path to be saved.
        """
        f, axes = plt.subplots(3, 2, figsize=size)
        plt.subplots_adjust(hspace=0.5)
        
        self.DrawCoV(axes[0, 0])
        self.DrawMbO(axes[1, 0])
        self.DrawInteraction(axes[2, 0])
        self.DrawMbP(axes[0, 1])
        self.DrawXbarChart(axes[1, 1])
        self.DrawRChart(axes[2, 1])
        
        if filename != None:
            plt.close()
            f.savefig(filename, dpi=200)
            
    def DrawRChart(self, ax):
        R = []
        for i in range(self.k):
            m = np.array([])
            temp = np.array(list(list(self.measurement.values())[i].values()))
            if len(m) == 0:
                m = np.max(temp,axis=1) - np.min(temp,axis=1)
            else:
                m = np.vstack((m, np.max(temp,axis=1) - np.min(temp,axis=1)))
            R.append(m)
        total = np.array(R).flatten()
        Rbar = np.mean(R)    
        
    #    UCL = Rbar + 3*np.std(R)
    #    LCL = max(Rbar - 3*np.std(R), 0)
        if self.k == 2:
            d2, d3 = 1.128, 0.8525
        elif self.k == 3:
            d2, d3 = 1.693, 0.8884
        elif self.k == 4:
            d2, d3 = 2.059, 0.8798
        elif self.k == 5:
            d2, d3 = 2.326, 0.8641
        elif self.k == 6:
            d2, d3 = 2.534, 0.8489
        elif self.k == 7:
            d2, d3 = 2.704, 0.8332
        elif self.k == 8:
            d2, d3 = 2.847, 0.8198
        elif self.k == 9:
            d2, d3 = 2.970, 0.8078
        elif self.k == 10:
            d2, d3 = 3.078, 0.7971
        else:
            raise ValueError('Not support the number of opperators.')
            
        UCL = Rbar + 3*d3*Rbar/d2
        LCL = max(Rbar - 3*d3*Rbar/d2, 0)   
        
        ax.set_title('R Chart')
        ax.set_xlabel('Operator')
        ax.set_ylabel('R')
        ax.hlines(xmin=0, xmax=len(total), y=Rbar, colors='g', linestyle='-')
        ax.text(len(total), Rbar, 'Rbar=%.2f' %(Rbar))
        ax.hlines(xmin=0, xmax=len(total), y=UCL, colors='r', linestyle='-')
        ax.text(len(total), UCL, 'UCL=%.2f' %(UCL))
        ax.hlines(xmin=0, xmax=len(total), y=LCL, colors='r', linestyle='-')
        ax.text(len(total), LCL, 'LCL=%.2f' %(LCL))
        
        r = np.max(total) - np.min(total)
        for i in range(len(R)):
            ax.vlines(x=i*len(R[i]), ymin=np.min(total)-0.1*r, ymax=np.max(total)+0.1*r, colors='b', linestyles='--')
        ax.set_xticks(np.arange(0, len(total), len(R[0])))
        ax.set_xticklabels(list(self.measurement.keys()))
        ax.plot(total, marker='o')
    
    def DrawXbarChart(self, ax):
        x = [] 
        data = []
        for i in range(len(self.measurement.keys())):
            m = np.array([])
            if len(m) == 0:
                m = np.array(list(list(self.measurement.values())[i].values()))
            else:
                m = np.vstack((m, np.array(list(list(self.measurement.values())[i].values()))))
            data.append(np.mean(m, axis=1))
            x.append(np.mean(m))
        Xbar = np.mean(x)
        UCL = np.mean(x) + 3*np.std(x)
        LCL = np.mean(x) - 3*np.std(x)
        total = np.array(data).flatten()
        
        ax.set_title('Xbar Chart')
        ax.set_xlabel('Operator')
        ax.set_ylabel('Xbar')    
        ax.hlines(xmin=0, xmax=len(total), y=Xbar, colors='g', linestyle='-')
        ax.text(len(total), Xbar, 'Xbar=%.2f' %(Xbar))
        ax.hlines(xmin=0, xmax=len(total), y=UCL, colors='r', linestyle='-')
        ax.text(len(total), UCL, 'UCL=%.2f' %(UCL))
        ax.hlines(xmin=0, xmax=len(total), y=LCL, colors='r', linestyle='-')
        ax.text(len(total), LCL, 'LCL=%.2f' %(LCL))
        
        r = np.max(total) - np.min(total)
        for i in range(len(data)):
            ax.vlines(x=i*len(data[i]), ymin=np.min(total)-0.1*r, ymax=np.max(total)+0.1*r, colors='b', linestyles='--')
        
        ax.set_xticks(np.arange(0, len(total), len(data[0])))
        ax.set_xticklabels(list(self.measurement.keys()))
        ax.plot(total, marker='o')  
    
    def DrawMbP(self, ax):
        ax.set_title('Measurement by Part')
        ax.set_xlabel('Part')
        ax.set_ylabel('Measurement')
        parts = list(self.measurement[list(self.measurement.keys())[0]].keys())
        part_dict = {}
        
        for i in parts:
            part_dict[i] = []    
        for app in self.measurement.keys():
            app_info = self.measurement[app]
            for part in app_info.keys():
                for p in app_info[part]:
                    part_dict[part].append(p)
        ax.boxplot(list(part_dict.values()))
        
        ax.set_xticklabels(parts)
        
    def DrawInteraction(self, ax):
        ax.set_title('Interaction Operator:Part')
        ax.set_xlabel('Part')
        ax.set_ylabel('mean of Measurement')
        
        x = list(self.measurement[list(self.measurement.keys())[0]].keys())
        ax.set_xticks(np.arange(min(x), max(x)+1, 1))
        
        data = []
        for i in range(len(self.measurement.keys())):
            m = np.array([])
            if len(m) == 0:
                m = np.array(list(list(self.measurement.values())[i].values()))
            else:
                m = np.vstack((m, np.array(list(list(self.measurement.values())[i].values()))))
            data.append(np.mean(m, axis=1).flatten())
        
        marker = itertools.cycle((',', '+', '.', 'o', '*')) 
        linestyle = itertools.cycle(('-', '--', '-.', ':', '*')) 
        for i in range(len(data)):        
            ax.plot(x, data[i], marker=next(marker), linestyle=next(linestyle), label='Op. %s' %(list(self.measurement.keys())[i]))
        ax.legend(loc=1)
        
    def DrawMbO(self, ax):
        ax.set_title('Measurement by Operator')
        ax.set_xlabel('Operator')
        ax.set_ylabel('Measurement')
        operators = self.measurement.keys()
        
        data = []
        for i in range(len(self.measurement.keys())):
            m = np.array([])
            if len(m) == 0:
                m = np.array(list(list(self.measurement.values())[i].values()))
            else:
                m = np.vstack((m, np.array(list(list(self.measurement.values())[i].values()))))
            data.append(m)
        ax.boxplot(data)    
        ax.set_xticklabels(operators)
        
    def DrawCoV(self, ax):
        ax.set_title('Components of Variation')
        xticks = ['totalRR', 'repeat.', 'reprod.', 'PartToPart']
        colors = ['red', 'lime', 'blue']
        
        width = 0.2
        ax.set_xticks(np.arange(len(xticks)) + width / 2)
        ax.set_xticklabels(xticks)
        
        if self.tolerance == 0 or self.tolerance == None:
            VarCompContrib = [self.GRRtable.iloc[0, 1],self.GRRtable.iloc[1, 1],self.GRRtable.iloc[2, 1],self.GRRtable.iloc[5, 1]]
            StudyCompContrib = [self.GRRtable.iloc[0, 4],self.GRRtable.iloc[1, 4],self.GRRtable.iloc[2, 4],self.GRRtable.iloc[5, 4]]
            
            rects0 = ax.bar(np.arange(len(xticks)) - 0.5*width, VarCompContrib, width, color=colors[0])
            rects1 = ax.bar(np.arange(len(xticks)) + 0.5* width, StudyCompContrib, width, color=colors[1])
            
            labels = ['VarCompContrib', 'StudyVarContrib']
            ax.legend((rects0[0], rects1[0]), labels, loc=1)            
        else:
            VarCompContrib = [self.GRRtable.iloc[0, 1],self.GRRtable.iloc[1, 1],self.GRRtable.iloc[2, 1],self.GRRtable.iloc[5, 1]]
            StudyCompContrib = [self.GRRtable.iloc[0, 4],self.GRRtable.iloc[1, 4],self.GRRtable.iloc[2, 4],self.GRRtable.iloc[5, 4]]
            ptRatio = [self.GRRtable.iloc[0, 5],self.GRRtable.iloc[1, 5],self.GRRtable.iloc[2, 5],self.GRRtable.iloc[5, 5]]
        
            rects0 = ax.bar(np.arange(len(xticks)) - 0.5*width, VarCompContrib, width, color=colors[0])
            rects1 = ax.bar(np.arange(len(xticks)) + 0.5* width, StudyCompContrib, width, color=colors[1])
            rects2 = ax.bar(np.arange(len(xticks)) + 1.5* width, ptRatio, width, color=colors[2])
           
            labels = ['VarCompContrib', 'StudyVarContrib', 'ptRatio']
            ax.legend((rects0[0], rects1[0], rects2[0]), labels, loc=1)