
class Chains():
    def __init__(self, iterations, ndim): #, lnprob, accept_rate):
        import numpy as np
        # initialize chain, acceptance rate and lnprob
        self.iterations = iterations
        self.ndim = ndim
        self.chain = np.zeros((iterations, ndim), dtype=np.float128)
        self.lnprob = np.zeros((iterations), dtype=np.float128)
        self.accept_rate = np.zeros((iterations), dtype=np.float128)
        #return 

    def __str__(self):
        return( "%s, %s, %s, %s, %s"
                 % (str(self.iterations), 
                    str(self.ndim),
                    str(self.chain),
                    #str(self.chain.shape),
                    str(self.lnprob),
                    str(self.accept_rate)) )
                    #str(self.lnprob.shape),
                    #str(self.accept_rate.shape)) )


    def dump(self, fname):
        import numpy as np
        fname_ch = "%s%s" % (fname, ".ch")
        fname_ar = "%s%s" % (fname, ".ar")
        fname_lp = "%s%s" % (fname, ".lp")
        np.savetxt(fname_lp, np.transpose([self.lnprob]))
        np.savetxt(fname_ar, np.transpose([self.accept_rate]))
        np.savetxt(fname_ch, self.chain)


