#!/usr/bin/env python
import os,glob,numpy
from wig import Wig
from rpy2.robjects.packages import importr
from rpy2.robjects import r,FloatVector
from copy import deepcopy
from time import time

'''
import sys,os,argparse,glob
'''
class Wigs:
    def __init__(self,path=None,step=0,suppress=False):
        '''
        Parameter:
            file: a pathe to the directory that contain the Wiggle format files
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region.
        '''
        self.data = {}  #in the format self.data[name]=Wig class instance
        self.step=step
        if path !=None:
            self.load(path=path,suppress=suppress)
    def ensureSameChrsByRemove(self):
        '''
        Description:
            make sure each Wig instance contain the same set of chrosomes, remove the chrosomes that are not contained by some instance
        Parameter:
            None
        Value:
            None
        '''
        wigs=self.data
        chrs={}
        for wig in wigs:
            for chr in wigs[wig].data:
               if chrs.has_key(chr):chrs[chr]+=1
               else:chrs[chr]=1
        crs=chrs.keys()
        wnum=len(wigs.keys())
        for cr in crs:
            if chrs[cr]<wnum:chrs.pop(cr)
        for wg in wigs:
            crs=wigs[wg].data.keys()
            for cr in crs:
                if not chrs.has_key(cr):wigs[wg].pop(cr)
        return wigs
    def foldNormalize(self,scalepairs=None,sampling_total=None,nonzero=False):
        '''
        Description:
            Normalize between Wig class instances by fold change
        Parameter:
            None
        Value:
            None
        '''
        ss=time()
        wigs=self.data
        names=wigs.keys()
        names.sort()
        
        if sampling_total==None:
            wsum={}
            for wig in wigs:wsum[wig]=wigs[wig].sum()
            asum=sum(wsum.values())/len(wsum.keys())
            for wig in names:
                print wig,'from',wigs[wig].sum(),'to',
                if scalepairs==None:wigs[wig].foldChange(asum*1.0/wsum[wig])
                else:wigs[wig].foldChange(scalepairs[wig]*1.0/wsum[wig])
                print wigs[wig].sum()
        else:
            average_total=sum(sampling_total.values())/len(sampling_total.keys())
            for name in names:
                print name,'from',wigs[name].sum(),'to',
                if scalepairs==None:wigs[name].foldChange(average_total*1.0/sampling_total[name])
                else:
                    wigs[name].foldChange(scalepairs[name]/sampling_total[name])
                print wigs[name].sum()
        if nonzero:
            print 'further correction based on count of non-zero base pairs'
            gsizes,non0sizes={},{}
            for wig in wigs:
                gsizes[wig]=wigs[wig].gsize()
                non0sizes[wig]=wigs[wig].non0size()
            agsize=sum(gsizes.values())*1.0/len(gsizes.keys())
            for wig in wigs:
                print wig,'from',wigs[wig].sum(),'to',
                wigs[wig].foldChange(non0sizes[wig]/agsize)
                print wigs[wig].sum(), 'based on non0size',non0sizes[wig],'and genome size',agsize
        return 1

    def samplingTotal(self,region_file=None,region_out_file=None,exclude_low_percent=1,exclude_high_percent=1,bnum=100000,nonzero=False):
        '''
        Description:
            caculate the sum of each wig's values after excluding the low and high percentile
        Parameter:
            None
        Value:
            None
        '''
        
        #if exclude_low_percent==0 and exclude_high_percent==0:return None
        #else:
        #print 'calculating normalization factors by sampling ...'
        names=self.data.keys()
        if exclude_low_percent==0 and exclude_high_percent==0 and region_file==None: return None
        sampling_total={}
        if region_file==None:
            print 'calculate total signal in each sample after excluding the top',exclude_high_percent,'and bottom',exclude_low_percent,'percents of genomic regions with extremely high and low signal values'
            wsums={}
            for name in names:wsums[name]=self.data[name].sum()
            wavg=sum(wsums.values())/len(wsums.values())
            
            rfwig=deepcopy(self.data[names[0]])
            rfwig.foldChange(wavg*1.0/wsums[names[0]])
            for name in names[1:]:
                self.data[name].foldChange(wavg*1.0/wsums[name])
                rfwig.add(self.data[name])
                self.data[name].foldChange(wsums[name]*1.0/wavg)
                
            rfwig.foldChange(1.0/len(names))
            lowcut,highcut=rfwig.percentile(p=[exclude_low_percent,100-exclude_high_percent],bnum=bnum,nonzero_end=nonzero)
            rg=rfwig.regionWithinValueRange(lowcut,highcut)
            if region_out_file!=None:rg.save(region_out_file)
        else:
            print 'calculate total signal in each sample in genomic regions defined by',region_file
            rg=Wig(region_file)
        for name in names:sampling_total[name]=self.data[name].multiply(rg).sum()
        print rg.sum(),'('+str(rg.sum()*100.0/rg.gsize())+'%) of',rg.gsize(),'base pairs calculated:'
        for name in names:print name,sampling_total[name],'('+str(sampling_total[name]*100.0/self.data[name].sum())+'% of total)'

        return sampling_total

    def get(self,k):
        '''
        Description:
            retrieve Wig class instance by name
        Parameter:
            k: the name of the Wig class instance
        Value:
            Wig class instance
        '''
        return self.data[k]
    def keys(self):
        '''
        Description:
            Retrieve the names of all Wig class instances
        Parameter:
            None
        Value:
            a list of names
        '''
        return self.data.keys()
    def load(self,path,suppress=False):
        '''
        Description:
            Load multiple Wig class instances from wiggle format files located in one directory
        Parameter:
            path: a path to the directory that contain the wiggle format files
        Value:
            None
        '''
        paths=path
        for path in paths.split(','):
            #wigs={}
            if os.path.isdir(path):
                for infile in glob.glob( os.path.join(path, '*.wig') ): 
                    fname=os.path.split(infile)[-1]
                    if fname[-4:]=='.wig':fname=fname[:-4]
                    self.set(fname,Wig(infile,step=self.step,suppress=suppress)) ########## ---add--- by kaifu on Aug 15,2012 ##########
                    #wigs[infile]=Wig(infile,step=self.step) ########## ---delete--- by kaifu on Aug 15,2012 ##########
            elif os.path.isfile(path):
                fname=os.path.split(path)[-1]
                if fname[-4:]=='.wig':fname=fname[:-4]
                self.set(fname,Wig(path,step=self.step,suppress=suppress)) ########## ---add--- by kaifu on Aug 15,2012 ##########
                #wigs[path]=Wig(path,step=self.step) ########## ---delete--- by kaifu on Aug 15,2012 ##########
            #self.data=wigs
    def nor(self,nor='F',exclude_low_percent=0,exclude_high_percent=0,scalepairs=None,sampling_total=None,nonzero=False):
        '''
        Description:
            normalize among Wig class instances.
        Parameter:
            nor: the normalization method, can be 'Q': quantile normalization, 'F':fold change scaling, 'S': sampling to same coverage, or 'N':no normalization
        Value:
            None
        '''
        if nor!='N' and len(self.data.keys())<2 and scalepairs==None:
            print "less than 2 datasets, no normalization to be done"
            return 1
        #if nor=='Q':
        #print "quantile normalization"
        #return self.quantileNormalize()
        elif nor=='F':
            #print "fold normalization"
            return self.foldNormalize(scalepairs=scalepairs,sampling_total=sampling_total,nonzero=nonzero)
        #elif nor=='S':
        #    return self.samplingNormalize()
        elif nor=='N':return 0
        #print("\nNo normalization method appointed to be done here.\n")
        else: print "Normalization method "+str(nor)+" not applicable now"
        return 0

    def pop(self,k):
        '''
        Description:
            remove Wig class instance by name
        Parameter:
            k: the name of the Wig class instance that is to be removed
        Value:
            None
        '''
        return self.data.pop(k)
    def quantileNormalize(self):
        '''
        Description:
            Normalize between Wig class instances by Quantile
        Parameter:
            None
        Value:
            None
        '''
        ss=time()
        self.ensureSameChrsByRemove()
        wigs=self.data
        r('require("preprocessCore")')
        normq=r('normalize.quantiles')
        
        chrs={}#now it is a dictionary, but will be change to a list later
        for wig in wigs:
            for chr in wigs[wig].data:
                if chrs.has_key(chr):chrs[chr]+=1
                else:chrs[chr]=1
        wnum=len(wigs.keys())
        '''
        pops=[]
        for chr in chrs:
            if chrs[chr]<wnum:pops.append(chr)
        for chr in pops:chrs.pop(chr)
        '''
        chrs=chrs.keys()#now chrs is a list
        names=wigs.keys()
        num=len(names)
        sizes={}
        size=0
        for chr in chrs:
            for name in names:
                if not wigs[name].data.has_key(chr):wigs[name].data[chr]=numpy.array([0.0])
                if not sizes.has_key(chr):sizes[chr]=wigs[name].data[chr].size
                elif sizes[chr]<wigs[name].data[chr].size:sizes[chr]=wigs[name].data[chr].size
            size+=sizes[chr]
        lst=numpy.array([0.0])
        lst.resize(size*num,refcheck=0)
        
        for i in range(0,num):
            name = names[i]
            tsize=size*i
            for j in range(0,len(chrs)):
                chr = chrs[j]
                wigs[name].data[chr].resize(sizes[chr],refcheck=0)
                ttsize=tsize+sizes[chr]
                lst[tsize:ttsize]+=wigs[name].data[chr][:sizes[chr]]
                tsize=ttsize
        mtr=r.matrix(FloatVector(lst),nrow = size, ncol = num)
        nmtr=normq(mtr)
        for i in range(0,num):
            name=names[i]
            tsize=size*i
            for j in range(0,len(chrs)):
                chr = chrs[j]
                ttsize=tsize+sizes[chr]
                wigs[name].data[chr][:sizes[chr]]=lst[tsize:ttsize]
                tsize=ttsize
        print 'time cost',time()-ss
        return 1

    def ajClonal(self,cut=1e-10,extend=1):  ###### add by Kaifu on Nov 14, 2012
        '''
        Description:
            Adjust clonal reads count, fold change between samples will not be altered in this process.
        Parameter:
            cut: the cutoff used to define clonal reads.
                When it is interger,  a read count larger than cut will be defined as clonal;
                when it is float, a read count that is larger than mean count by a Poisson test P value < cut will be defined as clonal.
            fsz: the extension length of each read that is used to calculate the wiggle file. Extension length means the length from 5' end to 3' end,
                e.g. a read may be 36bp when it is generated by the sequencing machine, but it might have been extended to be 80bp or cutted to be 1 bp,
                so the extension length will then be 80bp or 1bp.
        Value:
            None
        Note:
            all wiggle file in a Wigs object must have the same step size.
        '''
        
        #print '\nremoving clonal singal ...'
        ks=self.keys()
        m=deepcopy(self.get(ks[0]))
        if len(ks)>1:
            for k in ks[1:]:m.add(self.get(k))
        m.foldChange(1.0/len(ks))
        avg=m.mean()
        from math import log10,log
        if cut=='0':return
        else:
            if float(cut)>=1:cut=float(cut)
            else:
                co=cut.split('-')
                if len(co)==2:cut=float(co[1])-log10(float(co[0][:-1]))
                else:cut=0-log10(float(co[0]))
                lgpcut=cut
                cut=int(avg+0.5)
                ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log=TRUE))}''')
                while(0-(float(str(ppois(cut*1.0/extend,avg*1.0/extend)).split()[-1])/log(10))<lgpcut):cut+=1
        print 'aveage density is',avg,', use clonal signal cutoff',cut
        ks=self.keys()
        for chr in m.getChrs():
            tchrv=deepcopy(m.data[chr])
            tchrv-=cut#all positive values are count of clonal reads
            tchrv=((tchrv**2)**0.5+tchrv)/2#remove all neative values
            tchrv=m.data[chr]-tchrv+numpy.log(tchrv+1)# the addition of '1' is to avoid log(0),"+numpy.log(tchrv+1)" is used to keep the rank order values
            print chr+":"
            for k in ks:
                twg=self.get(k)
                if not twg.data.has_key(chr):continue
                temp=twg.data[chr].sum()
                print '\t'+k,'reduced from',temp,'to',
                if not twg.data.has_key(chr):twg.data[chr]=numpy.array([0.0])
                if twg.data[chr].size!=tchrv.size:twg.data[chr].resize(tchrv.size,refcheck=0)
                twg.data[chr]=tchrv*twg.data[chr]/(m.data[chr]+1e-100) # the addition of '1e-100' is to avoid devide by 0
                print twg.data[chr].sum(),', percent removed:',100-twg.data[chr].sum()*100.0/temp
        
        
    def samplingNormalize(self):
        '''
        Description:
            Normalize between Wig class instances by sampling to same coverage
        Parameter:
            None
        Value:
            None
        '''
        from random import randint
        ss=time()
        wigs=self.data
        wsum={}
        tarray=numpy.array([0.0])
        for wig in wigs:
            wsum[wig]=wigs[wig].sum()
        asum=sum(wsum.values())/len(wsum.keys())
        for k in wigs:
            wig=wigs[k]
            oldsum=wig.sum()
            num=oldsum-asum
            if num<0:num=0-num
            else:num=asum
            for chr in wig.data:
                tarray.resize(int(wig.data[chr].sum()),refcheck=0)
                tarray,csz,i,tsz=tarray*0,wig.data[chr].size,0,0
                while i<csz:
                    newtsz=int(tsz+wig.data[chr][i]+0.5)
                    if newtsz>=tarray.size:tarray.resize(newtsz+1000,refcheck=0)
                    while tsz<newtsz:
                        tarray[tsz]=i
                        tsz+=1
                    i+=1
                i,tnum,tsz=0,int((num*wig.chrSum(chr)/oldsum)/wig.step),tsz
                if oldsum>asum:wig.data[chr]*=0
                while i<tnum:
                    i+=1
                    wig.data[chr][tarray[randint(0,tsz-1)]]+=1
        print 'time cost',time()-ss
        return 1
    def set(self,k,wig):
        '''
        Description:
            add a Wig class instance
        Parameter:
            k: the name of the new Wig class instance
            wig: the new Wig class instance
        '''
        self.data[k]=wig
    def save(self,path,step=None,format="fixed"):
        '''
        Description:
            save all Wig class instances to wiggle format files in a directory.
            
        Parameter:
            path: the path to the directory that will contain the wiggle format files
            step: the stp size of the wiggle format files
            format: the format of the wiggle files, can be 'fixed' or 'var'
        Value:
            None
        '''
        if step==None:step=self.step
        wigs=self.data
        for k in wigs:
            if path!="":
                if not os.path.isdir(path):os.mkdir(path)
                df=os.path.split(str(k))
                tpath=os.path.join(path,df[-1])
            else:tpath=k
            if tpath[-4:]!='.wig':tpath=tpath+'.wig'
            wigs[k].save(file=tpath,format=format,step=step)
if __name__ == "__main__":
    print ''
    import sys

