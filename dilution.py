import numpy as np                 # Math library

class Dilution:
    
    def __init__(self, replica, strain, media, ref1, ref2):
        self.replica = replica
        self.strain = strain
        self.media = media
        self.ref1_bcid = ref1
        self.ref2_bcid = ref2
        self.ref1_idx = None
        self.ref2_idx = None
        self.ref1_counts = []
        self.ref2_counts = []
        self.BNSIZE = 5*(10**5)  # bottle neck size of dilution
        self.N_DILUTE = 8
        self.N_GEN = 5
        self.N_BCID = None
        self.bc_count_gen0 = []
        self.all_count = []
        self.pvals = [] # initial barcode frequecy
        self.Pvals = [] # sampling barcode freqency
        self.simu_reads = [] # simulated barcode reads
        #self.bcid = None
        #self.bcid_count = None
        #self.all_reads = None
        #self.record = {"BCID": [], "BCID_reads": [], "All_reads": []}
        
    def get_freq(self, bc_count, gen):
        return  bc_count/float(self.all_count[gen])
    
    def read_files(self):
        
        for gen in range(0,self.N_GEN):
            
            stitle = './source/' + self.replica + '_' + self.strain + '_' + self.media + '_gen' + repr(10*gen) + '.txt' # file name
            print stitle
            fopen = open(stitle,'r')  # open file
            dat = fopen.read().split()  # read file matrix and save in 1-dimension dat
            dat_rows = int(len(dat)/3)   # get length of coloum (number of row)
            dat = np.resize(dat,(dat_rows,3))  # resize the 1-dimension dat to matrix 
            fopen.close() # close file
            
            bc_count =[]
            bc_id = []
            
            for row in range(1,dat_rows):    # read rows, save first coloum as barcode index, second coloum as barcode counts
                bc_id.append(int(dat[row,0]))
                tmp = [map(float, x) for x in dat[row,1]] 
                bc_count.append(int(dat[row,1]))
                  
            bc_id = np.array(bc_id) # make format as numpy array
            bc_count = np.array(bc_count)  # make format as numpy array
            self.all_count.append(sum(bc_count)) # save all barcode counts
            
            itemindex = np.where(bc_id == self.ref1_bcid)  # find bc_count of the reference 
            self.ref1_idx = itemindex[0][0]
            self.ref1_counts.append(bc_count[itemindex[0][0]])
            itemindex = np.where(bc_id == self.ref2_bcid)
            self.ref2_idx = itemindex[0][0]
            self.ref2_counts.append(bc_count[itemindex[0][0]])
            if gen ==0:
                self.pvals = self.get_freq(bc_count,gen)
                self.N_BCID = len(self.pvals)
                #print bc_count
  
        print self.all_count
        print sum(self.pvals)
        
    def sampling(self, freq, n):
        sample = np.random.multinomial(n, freq, size=1)
        return sample[0]
    
    #def int_array(self, arr):
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
                self.simu_reads.extend(self.sampling(tmp_pvals,self.all_count[(dl-1)/2]))
               
        self.Pvals = np.resize(self.Pvals,(self.N_DILUTE,self.N_BCID))
        self.simu_reads = np.resize(self.simu_reads,(self.N_GEN,self.N_BCID)) 
        
        print self.Pvals[0]
        print (self.Pvals[0][self.ref1_idx], self.Pvals[0][self.ref2_idx]) # check    
            
#if __name__ == "__main__":
#    x = Dilution(REPLICA,STRAIN,MEDIA,BCID_REF1,BCID_REF2)