import pandas as pd
import pybedtools
from pybedtools.helpers import cleanup as pbtcleanup

class LoopRegion:
    def __init__(self, loopfn, regfn1, regfn2, loopid_col=None, reg1id_col=None, reg2id_col=None):
        reg1col = [0,1,2]
        if reg1id_col is not None:
            reg1col.append(reg1id_col)
        self.reg1 = pd.read_csv(regfn1, sep='\t', header=None)[reg1col]
        if self.reg1.shape[1]==3:
            self.reg1['regid'] = 'reg1_'+self.reg1.index.astype(str)
        self.reg1.columns = ['chr','start','end','regid']

        
        reg2col = [0,1,2]
        if reg2id_col is not None:
            reg2col.append(reg2id_col)
        self.reg2 = pd.read_csv(regfn2, sep='\t', header=None)[reg2col]
        if self.reg2.shape[1]==3:
            self.reg2['regid'] = 'reg2_'+self.reg2.index.astype(str)
        self.reg2.columns = ['chr','start','end','regid']        

        
        loopcol = [0,1,2,3,4,5]
        if loopid_col is not None:
            loopcol.append(loopid_col)
        self.loop = pd.read_csv(loopfn, sep='\t',header=None)[loopcol]
        if self.loop.shape[1]==6:
            self.loop['loopid'] = 'loop_'+self.loop.index.astype(str)
        self.loopl = self.loop[[0,1,2,'loopid']]
        self.loopl.columns = ['chr','start','end','loopid']
        self.loopr = self.loop[[3,4,5,'loopid']]
        self.loopr.columns = ['chr','start','end','loopid']

    
        reg1_bed = pybedtools.BedTool.from_dataframe(self.reg1)
        reg2_bed = pybedtools.BedTool.from_dataframe(self.reg2)
        loopl_bed = pybedtools.BedTool.from_dataframe(self.loopl)
        loopr_bed = pybedtools.BedTool.from_dataframe(self.loopr)
        
        
        self.int_l1 = loopl_bed.intersect(reg1_bed, wa=True, wb=True).to_dataframe()[['name','thickEnd']]
        self.int_l2 = loopl_bed.intersect(reg2_bed, wa=True, wb=True).to_dataframe()[['name','thickEnd']]
        self.int_r1 = loopr_bed.intersect(reg1_bed, wa=True, wb=True).to_dataframe()[['name','thickEnd']]
        self.int_r2 = loopr_bed.intersect(reg2_bed, wa=True, wb=True).to_dataframe()[['name','thickEnd']]
        self.int_l1.columns = ['loopid','reg1id']
        self.int_l2.columns = ['loopid','reg2id']
        self.int_r1.columns = ['loopid','reg1id']
        self.int_r2.columns = ['loopid','reg2id']
        
        pbtcleanup()
        
    def _reg1s_to_reg2s(self, reg1s):
        rtn = self.int_r2[
            self.int_r2['loopid'].isin(
                self.int_l1[
                    self.int_l1['reg1id'].isin(reg1s)
                ]['loopid']
            )
        ]['reg2id'].tolist() \
        + self.int_l2[
            self.int_l2['loopid'].isin(
                self.int_r1[
                    self.int_r1['reg1id'].isin(reg1s)
                ]['loopid']
            )
        ]['reg2id'].tolist()
        
        return list(set(rtn))
    
    def _reg2s_to_reg1s(self, reg2s):
        rtn = self.int_r1[
            self.int_r1['loopid'].isin(
                self.int_l2[
                    self.int_l2['reg2id'].isin(reg2s)
                ]['loopid']
            )
        ]['reg1id'].tolist() \
        + self.int_l1[
            self.int_l1['loopid'].isin(
                self.int_r2[
                    self.int_r2['reg2id'].isin(reg2s)
                ]['loopid']
            )
        ]['reg1id'].tolist()
    
        return list(set(rtn))

