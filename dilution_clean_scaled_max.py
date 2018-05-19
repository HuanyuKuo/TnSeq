import numpy as np                 # Math library
import matplotlib.pyplot as plt    # Plots
import sys                         # Directly print to output file
from scipy import stats            # linear regression


class Dilution:
    
    def __init__(self, replica, strain, media, ref, neutral, bnsize):
        self.BNSIZE = bnsize#1*(10**5)  # bottle neck size of dilution
        self.N_DILUTE = 8
        self.N_GEN = 5
        self.N_BCID = None
        #self.dat_rows = None
        self.S_REF = None
        self.replica = replica
        self.strain = strain
        self.media = media
        self.neutral = neutral
        self.ref_bcid =  ref
        self.ref_idx =  []# np.zeros((2,1))   
        self.ref_counts =  np.ones((2,self.N_GEN)) 
        self.pvals = [] # initial barcode frequecy (experimental value, and used as inital frequency in simualtion)
        self.Pvals = [] # sampling barcode freqency (simulated value)
        self.all_count = []  # all reads, experiment value
        self.simu_reads = [] # simulated barcode reads
        #self.exp_reads = None
        #self.exp_ini_freq = None
        self.sort_idx = [] # index of sorting of initial barcode freq (in decending order)
        self.s_record = None
        self.s_record_std = None
        #self.exp_exist = None
        #self.s_dat = []
        self.s_bcid = []
        self.s_inifreq = []
        self.s_endfreq =[]
        self.s_freq = []
        self.s_slope = []
        self.s_intercept = []
        self.s_rvalue = []
        self.s_pvalue = []
        self.s_stderr = []
        self.s_var =[]
        #self.s_r2 = []
        self.tmp=None
        self.bc_count = None
        self.bc_id = []
        
    def runs(self):
        # import inital bc_freq from file
        #self.read_files()     
        # simulate dilution and reads
        self.dilutions()
        # calculate simulated selection coefficient
        simu_scaled_reads = self.rescale_reads(self.simu_reads, self.ref_idx[0], self.ref_idx[1])
        self.s_record, self.s_record_std = self.selection_coef(simu_scaled_reads)
    
    def get_experimental_s(self):
        self.s_bcid = []
        self.s_inifreq = []
        self.s_slope = []
        self.s_intercept = []
        self.s_rvalue = []
        self.s_pvalue = []
        self.s_stderr = []
        for bc in range(0,len(self.bc_id)):
            y = []
            x = []
            for gen in range(0,self.N_GEN):
                if self.bc_count[gen][bc]!= 0:
                    y.append(float(self.bc_count[gen][bc])/float(self.ref_counts[0][gen]+self.ref_counts[1][gen]))#float(self.all_count[gen]))
                    x.append(10.0*gen)
            y = np.log(y)
            #print self.bc_id[bc],x,y
            if len(x) >=5:
                self.s_bcid.append(self.bc_id[bc])
                self.s_inifreq.append(self.pvals[bc])
                self.s_endfreq.append(float(self.bc_count[4][bc])/float(self.all_count[4]))
                tmp = [float(self.bc_count[gen][bc])/float(self.all_count[gen]) for gen in range(0,5) ]
                self.s_freq.append(tmp)
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                self.s_slope.append(slope)
                self.s_intercept.append(intercept)
                self.s_rvalue.append(r_value)
                self.s_pvalue.append(p_value)
                self.s_stderr.append(std_err)
                self.s_var.append(np.var(y))

       
    def get_index(self,arr,val):
        index = np.flatnonzero(arr == val)
        if not index: # if is empty
            return float('NaN')
        else:
            return index[0]
    
    def read_files(self):
       
        # Read input file:
        filename = './source/' + self.replica + '_' + self.strain + '_' + self.media + '_clean.txt' # file name
        fopen = open(filename,'r')  # open file
        dat = fopen.read().split('\n')  # read file matrix and save in 1-dimension dat
        fopen.close() # close file
        dat2 = []
        COL = 10
        ROW = len(dat)
        for row in range(0, ROW):
            tmp = dat[row].split('\t')
            dat2.extend(tmp)
        dat2 = np.resize(dat2,(ROW,COL))  # resize the 1-dimension dat to matrix 
        dat = dat2[1:-1]
        ROW = len(dat)
        del dat2
        # Read barcode index and barcode counts:
        bc_count =[]
        bc_id = []
        bc_id.extend( [int(dat[row][0]) for row in range(0,ROW)])
        bc_count.extend( [  [int(dat[row][gen+5]) for row in range(0,ROW) ]]for gen in range(0,self.N_GEN))
        tmp = [bc_count[gen][0] for gen in range(0,self.N_GEN)]
        bc_count = tmp
        del tmp
        #for row in range(0,ROW):
           # if dat[row][1] == '1019':
           #     print dat[row]
        self.tmp = dat
        bc_id = np.array(bc_id) # make format as numpy array
        bc_count = np.array(bc_count)  # make format as numpy array
        self.all_count = ([sum(bc_count[gen]) for gen in range(0,self.N_GEN)]) # save all barcode counts
        # Record initial barcode freqnecy as "real" probability for samping, make it in decending order
        pvals = bc_count[0]/float(sum(bc_count[0]))
        self.sort_index = np.argsort(bc_count[0]) 
        self.sort_index = self.sort_index[::-1] # reverse
        self.pvals = pvals[self.sort_index] # make it in decending order
        self.bc_id = bc_id[self.sort_index]
        self.bc_count = [bc_count[gen][self.sort_index] for gen in range(0, self.N_GEN)]
        self.N_BCID = len(self.pvals)
        # Find bc_count & idx of the reference:
        for ref in range(0,2):
                idx1 = self.get_index(bc_id, self.ref_bcid[ref])  # index of reference barcode in original file
                idx2 = self.get_index(self.sort_index, idx1) # index of reference barcode in sorted indxing
                self.ref_idx.append(idx2)
                for gen in range(0, self.N_GEN):
                    #print idx1, idx2
                    self.ref_counts[ref][gen] = self.bc_count[gen][idx2]
                    
 
    def sampling(self, freq, n):
        sample = np.random.multinomial(n, freq, size=1)
        return sample[0]
    
    def f(self,x):
        return np.int(x)
    
    def dilutions(self):
        
        # before generation zero:    
        get_int = np.vectorize(self.f) 
        #self.Pvals.extend(self.pvals)
        self.simu_reads.extend(get_int(self.pvals*self.all_count[0]))
        
        # start diluting and growing:
        tmp_pvals = self.pvals
        for dl in range(0,self.N_DILUTE): # 8 dilutions cause sampling fluctuation
            tmp_pvals = self.sampling(tmp_pvals,self.BNSIZE)/float(self.BNSIZE)
            self.Pvals.extend(tmp_pvals)
            if dl%2 == 1:  # 4 reads cause sampling fluctuation
                self.simu_reads.extend(self.sampling(tmp_pvals,self.all_count[(dl+1)/2]))
               
        self.Pvals = np.resize(self.Pvals,(self.N_DILUTE,self.N_BCID)) # deceding order as pvals
        self.simu_reads = np.resize(self.simu_reads,(self.N_GEN,self.N_BCID)) # decending order as pvals
        
        # Force reference reads as experimental value:
        for ref in range(0,2):
            idx = self.ref_idx[ref]
            for gen in range(0,self.N_GEN):
                self.simu_reads[gen,idx] = self.ref_counts[ref][gen]
        # Get selection coeff. of refernece strain:
        self.get_s_REF()
                
    def rescale_reads(self,reads,idx1,idx2):
        scaled_reads =[]
        for gen in range(0,self.N_GEN):
            reference_reads = reads[gen,idx1]+reads[gen,idx2]
            scaled_reads.extend(reads[gen,0:self.N_BCID] / float(reference_reads))
            
        scaled_reads = np.resize(scaled_reads,(self.N_GEN,self.N_BCID))
        return scaled_reads
    
    def rescale_reads_exp_max(self,reads,idx1,idx2):
        scaled_reads =[]
        for gen in range(0,self.N_GEN):
            #reference_reads = reads[gen,idx1]+reads[gen,idx2]
            reference_reads = reads[gen][0] # max ini freq
            #reference_reads = self.bc_count[gen][0] # max ini freq
            scaled_reads.extend(reads[gen,0:self.N_BCID] / float(reference_reads))
            
        scaled_reads = np.resize(scaled_reads,(self.N_GEN,self.N_BCID))
        return scaled_reads
    def selection_coef(self, scaled_reads):
        
        s_simu = [0 for i in range(0,self.N_BCID)]
        s_simu_std = [0 for i in range(0,self.N_BCID)]
        x = np.arange(0,45,10)
        for bc in range(0,self.N_BCID):
            y = np.log(scaled_reads[0:self.N_GEN,bc])
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            
            #z = np.polyfit(x, y, 1)
            
            s_simu[bc] = slope # z[0]
            s_simu_std[bc] = std_err
        return s_simu, s_simu_std
        
    # -------------------#
    # Optional Functions #
    # -------------------#
   
    def output_simu_reads(self): # Optional
        
        stitle = './source/' + self.replica + '_' + self.strain + '_' + self.media + '_gen0' + '.txt' # file name
        
        dat, rows = self.input_dat(stitle,3)
        
        stitle = './output/' + self.replica + '_' + self.strain + '_' + self.media + '_simulated_counts' + '.csv' # output file name
        fout = open(stitle,'w')
        orig_stdout = sys.stdout  #using system to 'print' to ouput file
        sys.stdout = fout
        #print dat
        print 'STRAIN\tBC_ID\tMU_ID\tENV\tgen_0\tgen_10\tgen_20\tgen_30\tgen_40'
        for idx_BCID in range(0,self.N_BCID):
            outstring =   self.strain + '\t' + dat[1+idx_BCID,0] + '\t' + dat[1+idx_BCID,2] + '\t' + self.replica + '_' + self.media
            for idx_gen in range(0,5):
                bcid = int(dat[1+idx_BCID,0]) # get value of bacarod index 
                idx_sorted = self.get_index(self.sort_bcid0, bcid) # search the sorted index for this bcid
                outstring += '\t' + repr(self.simu_reads[idx_gen,idx_sorted])
            print outstring
        sys.stdout = orig_stdout
        fout.close()
    
    def sort_s_for_plot(self, s_record):
        
        sort_indx = np.argsort(s_record)
        #sort_indx = sort_indx[::-1]  
        tmp = np.array(s_record)
        tmp = tmp[sort_indx]
        itemindex = np.flatnonzero(np.isnan(tmp))
        if len(itemindex)==0: # if it's empty
            UPBOUND = len(tmp)
        else:
            UPBOUND = itemindex[0]
        s_plot = tmp[0:UPBOUND]
        return UPBOUND, s_plot
       
    def plot_s_dist(self, s_record):  # Optional
        
        UPBOUND, s_plot = self.sort_s_for_plot(s_record)
        print UPBOUND
        width = 2*np.sqrt(np.var(s_plot))
        mean = np.mean(s_plot)
        sN = "N {0:0d}".format(len(s_plot))
        sw = "Width {0:.4f}".format(width)
        sm = "Mean {0:0.5f}".format(mean)
        #sr = "Reference strain slope {0:0.5f}".format(self.S_REF)
        plt.figure()
        s_percent = [100*s_plot[i] for i in range(0,UPBOUND)]
        slabel = self.strain + '_' + self.replica + '_' + self.media
        plt.hist(s_percent, 16, label = slabel)
        plt.title('Histogram of selection coefficient')
        plt.xlabel('s (per generation) %')
        plt.legend(loc="best")
        plt.show()
        print (sN+'\n'+sw+'\n'+sm+'\n')#+sr)
       
    def plot_scale_freq(self): # Optional
        plt.figure(num=None, figsize=(15, 5), dpi=80)
        plt.plot([2,1,3])
        plt.subplot(121)
        x = np.arange(0,45,10)
        col = plt.cm.jet(np.linspace(0,1,20)) 
        scaled_reads = self.rescale_reads(self.simu_reads, self.ref_idx[0], self.ref_idx[1])
        for bc in range(0,20):#self.N_BCID/100):
            s = self.s_record[bc] 
            if np.isnan(s) == False: # if s is NOT Nan
                y = np.log(scaled_reads[0:self.N_GEN,bc])
                cidx = int( (((0.1+s)>0)&((-0.1+s)<0))*(0.1+s)*100+ ((-0.1+s)>=0)*10)
                #plt.plot(x,y,color=col[cidx],marker='.')
                plt.plot(x,y,color=col[bc],marker='.')
        plt.title("SIMU Log of BC Freqeucy scaled by reference strains")
        plt.xlabel('generation')
        plt.ylim(0.6,2.5)
        plt.xlim(-1,41)
        #plt.show()
        
        plt.subplot(122)
        x = np.arange(0,45,10)
        col = plt.cm.jet(np.linspace(0,1,20)) 
        #scaled_reads = self.rescale_reads(self.simu_reads, self.ref_idx[0], self.ref_idx[1])
        #self.s_freq = []

        for bc in range(0,20):#self.N_BCID/100):
            s = self.s_slope[bc] 
            if np.isnan(s) == False: # if s is NOT Nan
                tmp = np.multiply(self.s_freq[bc][0:self.N_GEN],self.all_count)
                ref_count = np.add(self.ref_counts[0][0:self.N_GEN],self.ref_counts[1][0:self.N_GEN])
                tmp = np.divide(tmp,ref_count)
                y = np.log(tmp)
                cidx = int( (((0.1+s)>0)&((-0.1+s)<0))*(0.1+s)*100+ ((-0.1+s)>=0)*10)
                plt.plot(x,y,color=col[bc],marker='.')
        plt.title("EXP Log of BC Freqeucy scaled by reference strains")
        plt.xlabel('generation')
        plt.ylim(0.6,2.5)
        plt.xlim(-1,41)
        plt.show()
    def get_s_REF(self):        
        x2 = np.arange(0,45,10)
        y2_all = [0 for i in range(0,self.N_GEN)]
        for ref in range(0,2):
            for gen in range(0,self.N_GEN):
                idx = self.ref_idx[ref]
                bc_freq = float(self.simu_reads[gen,idx])/float(self.all_count[gen])
                y2_all[gen] += bc_freq       
        y2_all = np.log(y2_all)
        z = np.polyfit(x2, y2_all, 1)
        fitline = np.poly1d(z)
        self.S_REF = z[0]
        
            
    def plot_ref_strain(self):  # Optional
       
        # Plot two figures of reference- barcode freqency:
        x = np.arange(5,45,5)
        x2 = np.arange(0,45,10)
        x3 = np.arange(0,45,1)
        y2 =[0 for i in range(0,self.N_GEN)]
        y_all = [0 for i in range(0,self.N_DILUTE)]
        y2_all = [0 for i in range(0,self.N_GEN)]
        for ref in range(0,2):
            plt.figure()
            # Get (x,y) for dilution, (x2,y2) for reads, (x3, y3) for fit line
            bc_freq = self.Pvals[0:self.N_DILUTE,self.ref_idx[ref]]
            y = np.log(bc_freq)        
            y_all += bc_freq
            for gen in range(0,self.N_GEN):
                idx = self.ref_idx[ref]
                bc_freq = float(self.simu_reads[gen,idx])/float(self.all_count[gen])
                y2[gen] = np.log(bc_freq)
                y2_all[gen] += bc_freq
            # Get fit line
            z = np.polyfit(x2, y2, 1)
            fitline = np.poly1d(z)
            y3 = fitline(x3)
            # make plots
            plt.plot(x,y,'bo',label='Dilution')
            plt.plot(x2,y2,'ro',label='Reads')
            plt.plot(x3,y3,'r',label='SLOPE {0:.3f}'.format(z[0]))
            plt.xlabel('generation')
            plt.legend()
            plt.title('Log of BC freqency of Referece strain '+ self.neutral +' {0:0d}'.format(self.ref_bcid[ref]))
            plt.show()
            
        # Plot reference together as one figure:
        plt.figure()
        # Get y
        y_all = np.log(y_all)
        y2_all = np.log(y2_all)
        # Get fit line
        z = np.polyfit(x2, y2_all, 1)
        fitline = np.poly1d(z)
        y3_all = fitline(x3)
        #self.S_REF = z[0]
        # make plot
        plt.plot(x,y_all,'bo',label='Dilution')
        plt.plot(x2,y2_all,'ro',label='Reads')
        plt.plot(x3,y3_all,'r',label='SLOPE {0:.3f}'.format(self.S_REF))
        plt.xlabel('generation')
        plt.legend()
        plt.title('Log of BC freqency of both Referece strains.')
        plt.show()
       
   
        
       

   
        
          
     