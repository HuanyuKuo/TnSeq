import numpy as np                 # Math library
import matplotlib.pyplot as plt    # Plots
import sys                         # Directly print to output file


class Dilution:
    
    def __init__(self, replica, strain, media, ref, neutral):
        self.BNSIZE = 5*(10**5)  # bottle neck size of dilution
        self.N_DILUTE = 8
        self.N_GEN = 5
        self.N_BCID = None
        self.dat_rows = None
        self.S_REF = None
        self.replica = replica
        self.strain = strain
        self.media = media
        self.neutral = neutral
        self.ref_bcid =  ref
        self.ref_idx =  np.ones((2,self.N_GEN))   
        self.ref_counts =  np.ones((2,self.N_GEN)) 
        self.pvals = [] # initial barcode frequecy
        self.Pvals = [] # sampling barcode freqency
        self.all_count = []  # all reads, experiment value
        self.simu_reads = [] # simulated barcode reads
        self.scaled_reads = [] # rescale simlated barcode reads to reference value
        self.sort_reads = []
        self.s_record = None
        
    def runs(self):
        self.read_files()
        self.dilutions()
        self.selection_coef()
        # Optional Outputs and Display
        #self.output_simu_reads()
        #self.plot_ref_strain()  
       
    def get_freq(self, bc_count, gen):
        return  bc_count/float(self.all_count[gen])
    
    def input_dat(self, file):
        fopen = open(file,'r')  # open file
        dat = fopen.read().split()  # read file matrix and save in 1-dimension dat
        self.dat_rows = int(len(dat)/3)   # get length of coloum (number of row)
        dat = np.resize(dat,(self.dat_rows,3))  # resize the 1-dimension dat to matrix 
        fopen.close() # close file
        return dat
    
    def get_index(self,arr,val):
        itemindex = np.where(arr == val)
        return itemindex[0][0]
    
    def read_files(self):
        
        for gen in range(0,self.N_GEN):
            
            # Read input file:
            filename = './source/' + self.replica + '_' + self.strain + '_' + self.media + '_gen' + repr(10*gen) + '.txt' # file name
            dat = self.input_dat(filename)
            
            # Read rows. Save first & second coloum as barcode index and barcode counts:
            bc_count =[]
            bc_id = []
            
            for row in range(1,self.dat_rows):   
                bc_id.append(int(dat[row,0]))
                tmp = [map(float, x) for x in dat[row,1]] 
                bc_count.append(int(dat[row,1]))
                  
            bc_id = np.array(bc_id) # make format as numpy array
            bc_count = np.array(bc_count)  # make format as numpy array
            self.all_count.append(sum(bc_count)) # save all barcode counts
            
            # Find bc_count & idx of the reference:
            for ref in range(0,1):
                idx = self.get_index(bc_id, self.ref_bcid[ref])  
                self.ref_idx[ref][gen]= idx
                self.ref_counts[ref][gen] = bc_count[idx]
            
            # Record initial barcode freqnecy as "real" probability for samping
            if gen ==0:
                self.pvals = self.get_freq(bc_count,gen)
                self.N_BCID = len(self.pvals)
        
    def sampling(self, freq, n):
        sample = np.random.multinomial(n, freq, size=1)
        return sample[0]
    
    def f(self,x):
        return np.int(x)
    
    def dilutions(self):
        
        # before generation zero:    
        f2 = np.vectorize(self.f)
        self.Pvals.extend(self.pvals)
        self.simu_reads.extend(f2(self.pvals*self.all_count[0]))
        
        # start diluting and growing:
        tmp_pvals = self.pvals
        for dl in range(0,self.N_DILUTE): # 8 dilutions cause sampling fluctuation
            tmp_pvals = self.sampling(tmp_pvals,self.BNSIZE)/float(self.BNSIZE)
            self.Pvals.extend(tmp_pvals)
            if dl%2 == 1:  # 4 reads cause sampling fluctuation
                self.simu_reads.extend(self.sampling(tmp_pvals,self.all_count[(dl+1)/2]))
               
        self.Pvals = np.resize(self.Pvals,(self.N_DILUTE,self.N_BCID))
        self.simu_reads = np.resize(self.simu_reads,(self.N_GEN,self.N_BCID))
        
    def rescale_reads(self):
        ref1 = self.ref_idx[0][0]
        ref2 = self.ref_idx[1][0]
        for gen in range(0,self.N_GEN):
            reference_reads = self.simu_reads[gen,ref1]+self.simu_reads[gen,ref2]
            self.scaled_reads.extend(self.simu_reads[gen,0:self.N_BCID] / float(reference_reads))
            
        self.scaled_reads = np.resize(self.scaled_reads,(self.N_GEN,self.N_BCID))
        
    def sort_scaled_reads(self):
        
        self.rescale_reads()
        
        logcheck = np.sum(self.scaled_reads,axis=0)
        logcheck_sort_indx = np.argsort(logcheck)
        logcheck_sort_indx = logcheck_sort_indx[::-1]  # reverse
        self.sort_reads = self.scaled_reads[0:self.N_GEN,logcheck_sort_indx]
        itemindex = np.where(self.sort_reads[self.N_GEN-1,0:self.N_BCID]==0)
        self.UPBOUND = itemindex[0][0]
        self.sort_reads = self.sort_reads[0:self.N_GEN,0:self.UPBOUND]
       
        
    def selection_coef(self):
        self.sort_scaled_reads()
        self.s_record = [0 for i in range(0,self.UPBOUND)]
        x = np.arange(0,45,10)
        for bc in range(0,self.UPBOUND):
            y = np.log(self.sort_reads[0:self.N_GEN,bc])
            z = np.polyfit(x, y, 1)
            self.s_record[bc] = z[0]
        
    # -------------------#
    # Optional Functions #
    # -------------------#
   
    def output_simu_reads(self): # Optional
        
        stitle = './source/' + self.replica + '_' + self.strain + '_' + self.media + '_gen0' + '.txt' # file name
        dat = self.input_dat(stitle)
        
        stitle = './output/' + self.replica + '_' + self.strain + '_' + self.media + '_simulated_counts' + '.csv' # output file name
        fout = open(stitle,'w')
        orig_stdout = sys.stdout  #using system to 'print' to ouput file
        sys.stdout = fout
        #print dat
        print 'STRAIN\tBC_ID\tMU_ID\tENV\tgen_0\tgen_10\tgen_20\tgen_30\tgen_40'
        for idx_BCID in range(0,self.N_BCID):
            outstring =   self.strain + '\t' + dat[1+idx_BCID,0] + '\t' + dat[1+idx_BCID,2] + '\t' + self.replica + '_' + self.media
            for idx_gen in range(0,5):
                outstring += '\t' + repr(self.simu_reads[idx_gen,idx_BCID])
            print outstring
        sys.stdout = orig_stdout
        fout.close()
        
    def plot_s_dist(self):  # Optional
        width = 2*np.sqrt(np.var(self.s_record))
        mean = np.mean(self.s_record)
        sN = "N {0:0d}".format(len(self.s_record))
        sw = "Width {0:.4f}".format(width)
        sm = "Mean {0:0.5f}".format(mean)
        sr = "Reference strain slope {0:0.5f}".format(self.S_REF)
        plt.figure()
        s_percent = [100*self.s_record[i] for i in range(0,self.UPBOUND)]
        slabel = self.strain + '_' + self.replica + '_' + self.media
        plt.hist(s_percent, 16, label = slabel)
        plt.title('Histogram of selection coefficient')
        plt.xlabel('selection coefficient (per generation) %')
        plt.legend()
        plt.show()
        print (sN+'\n'+sw+'\n'+sm+'\n'+sr)
       
    def plot_scale_freq(self): # Optional
        plt.figure()
        x = np.arange(0,45,10)
        col = plt.cm.jet(np.linspace(0,1,20)) 
        for bc in range(0,self.UPBOUND):
            y = np.log(self.sort_reads[0:self.N_GEN,bc])
            s = self.s_record[bc] 
            cidx = int( (((0.1+s)>0)&((-0.1+s)<0))*(0.1+s)*100+ ((-0.1+s)>=0)*10)
            plt.plot(x,y,color=col[cidx])
        plt.title("Log of BC Freqeucy scaled by reference strains")
        plt.xlabel('generation')
        plt.show()
                   
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
            bc_freq = self.Pvals[0:self.N_DILUTE,self.ref_idx[ref][0]]
            y = np.log(bc_freq)        
            y_all += bc_freq
            for gen in range(0,self.N_GEN):
                idx = self.ref_idx[ref][0]
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
        self.S_REF = z[0]
        # make plot
        plt.plot(x,y_all,'bo',label='Dilution')
        plt.plot(x2,y2_all,'ro',label='Reads')
        plt.plot(x3,y3_all,'r',label='SLOPE {0:.3f}'.format(self.S_REF))
        plt.xlabel('generation')
        plt.legend()
        plt.title('Log of BC freqency of both Referece strains.')
        plt.show()
       
   
        
       

   
        
          
     