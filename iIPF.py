# seed, colsum, rowsum
# should be sum(colsum) = sum(rowsum)
#

class IPF_2dInteger(object):

    def __init__(self, seed, rowsum, colsum):

        if abs(sum(colsum) - sum(rowsum)) > 0.1:
            return None

        # self.seed = seed
        self.colsum = colsum
        self.rowsum = rowsum
        self.ncol   = len(colsum)
        self.nrow   = len(rowsum)
        # self.total  = sum(rowsum)

        self.A = seed[:][:]
        self.intA = []
        self.fA = []
        self.report = False
        # self.icolsum = []
        # self.irowsum = []

        # getting ready the seed for IPF : clarify the meaning of zero.
        for i in xrange(self.nrow):
            for j in xrange(self.ncol):
                if self.A[i][j] == 0:
                    self.A[i][j] = 1.0e-10

        for i in xrange(self.nrow):
            if rowsum[i] == 0:
                for j in xrange(self.ncol):
                    self.A[i][j] = 0.0

        for j in xrange(self.ncol):
            if colsum[j] == 0:
                for i in xrange(self.nrow):
                    self.A[i][j] = 0.0
        return


    def fIPF0(self):
        return self.fIPF(self.A, self.rowsum, self.colsum)

    # generic ipf
    def fIPF(self, matrix, frow, fcol):

        for iter in xrange(100):

            maxfactor_r = -1.0e10
            minfactor_r =  1.0e10
            for i in xrange(self.nrow):
                if frow[i] == 0:
                    continue

                sumx = 0.0
                for j in xrange(self.ncol):
                    sumx += matrix[i][j]

                factor = frow[i] / sumx
                if maxfactor_r < factor: maxfactor_r = factor
                if minfactor_r > factor: minfactor_r = factor

                for j in xrange(self.ncol):
                    matrix[i][j] *= factor

            maxfactor_c = -1.0e10
            minfactor_c =  1.0e10
            for j in xrange(self.ncol):
                if fcol[j] == 0:
                    continue

                sumy = 0.0
                for i in xrange(self.nrow):
                    sumy += matrix[i][j]

                factor = fcol[j] / sumy
                if maxfactor_c < factor: maxfactor_c = factor
                if minfactor_c > factor: minfactor_c = factor

                for i in xrange(self.nrow):
                    matrix[i][j] *= factor

            if abs(maxfactor_r - 1.0) < 1.0e-8 and abs(1.0 - minfactor_r) < 1.0e-8 and \
                abs(maxfactor_c - 1.0) < 1.0e-8 and abs(1.0 - minfactor_c) < 1.0e-8:
                break

        return iter


    @staticmethod
    def real_to_integer_and_float(matrix):
        imat = []
        fmat = []

        for rowvector in matrix:
            itmp, ftmp = [], []
            for value in rowvector:
                ivalue = int(value)
                itmp.append(ivalue)
                ftmp.append(value - ivalue)
            imat.append(itmp)
            fmat.append(ftmp)

        return imat, fmat

    # the sum always should be integer
    @staticmethod
    def marginals(matrix):
        nrow = len(matrix)
        ncol = len(matrix[0])

        sumrow = [0] * nrow
        sumcol = [0] * ncol
        for i in xrange(nrow):
            for j in xrange(ncol):
                sumrow[i] += matrix[i][j]
                sumcol[j] += matrix[i][j]


        for i in xrange(nrow):
            sumrow[i] = int(sumrow[i] + 0.1)

        for j in xrange(ncol):
            sumcol[j] = int(sumcol[j] + 0.1)
        return sumrow, sumcol


    def nonezerocells_to_ones(self, frow, fcol):
        # check if the number of nonzero cells is same to the remaining row sum...

        tmpMatrix = []
        for i in xrange(self.nrow):
            tmpMatrix.append([0]*self.ncol)
            nzero = 0
            for j in xrange(self.ncol):
                if self.fA[i][j] > 0:
                    nzero += 1

            if abs(1.0 * nzero - frow[i]) < 1.0e-8:
                for j in xrange(self.ncol):
                    if self.fA[i][j] > 0:
                        tmpMatrix[i][j] = 1

        for j in xrange(self.ncol):
            nzero = 0
            for i in xrange(self.nrow):
                if self.fA[i][j] > 0:
                    nzero += 1

            if abs(1.0 * nzero - fcol[j]) < 1.0e-8:
                for i in xrange(self.nrow):
                    if self.fA[i][j] > 0:
                        tmpMatrix[i][j] = 1

        sumfA = 0
        for i in xrange(self.nrow):
            for j in xrange(self.ncol):
                if tmpMatrix[i][j] == 1:
                    self.intA[i][j] += 1
                    self.fA[i][j] = 0.0
                    frow[i] -= 1
                    fcol[j] -= 1
                sumfA += self.fA[i][j]

        sumfrow = 0
        for i in xrange(self.nrow):
            sumfrow += frow[i]

        sumfcol = 0
        for j in xrange(self.ncol):
            sumfcol += fcol[j]

        return sumfA, sumfrow, sumfcol


    def drop_minvalues(self):
        minvalue = 1.0e10
        matrix = self.fA
        for i in xrange(self.nrow):
            for j in xrange(self.ncol):
                if matrix[i][j] > 0:
                    if minvalue > matrix[i][j]:
                        minvalue = matrix[i][j]

        if minvalue < 0.01: minvalue = 0.01

        # check whether the matrix would have anything after setting the min values to zero
        sumNotMin = 0
        for i in xrange(self.nrow):
            for j in xrange(self.ncol):
                if matrix[i][j] > minvalue:
                    sumNotMin += matrix[i][j]
        if sumNotMin == 0:
            return

        for i in xrange(self.nrow):
            for j in xrange(self.ncol):
                if matrix[i][j] <= minvalue:
                    matrix[i][j] = 0
        return


    def rollup_maxvalues(self, frow, fcol):
        matrix = self.fA
        xsum = 0
        for i in xrange(self.nrow):
            for j in xrange(self.ncol):
                if matrix[i][j] >= 1.00:
                    if self.report:
                        if matrix[i][j] > 1.99:
                            print 'too concentrated %d %d %f' % (i, j, matrix[i][j])

                    if frow[i] >= 1.0 and fcol[j] >= 1.0:
                        matrix[i][j] = 0
                        self.intA[i][j] += 1
                        frow[i] -= 1
                        fcol[j] -= 1
                        xsum += 1
        return xsum

    # a method that matches the colsums by the bucket roundings
    def iIPF_column_base(self):
        main_iteration = self.fIPF(self.A, self.rowsum, self.colsum)

        imat = []
        for i in range(self.nrow):
            imat.append([0] * self.ncol)

        for j in range(self.ncol):
            if self.colsum[j] == 0: continue
            sumfcol = 0.0
            for i in range(self.nrow):
                sumfcol += self.A[i][j]

            rate = self.colsum[j] / sumfcol
            residual = 0.0
            for i in range(self.nrow):
                fvalue = self.A[i][j] * rate
                ivalue = int(fvalue)
                if residual + fvalue - ivalue >= 0.5:
                    ivalue += 1
                residual += fvalue - ivalue
                imat[i][j] = ivalue

        self.intA = imat
        return

    def iIPF_row_base(self):
        main_iteration = self.fIPF(self.A, self.rowsum, self.colsum)
        
        imat = []
        for i in range(self.nrow):
            imat.append([0] * self.nrow)
            
            if self.rowsum[i] == 0: continue
            sumfrow = 0.0
            for j in range(self.ncol):
                sumfrow += self.A[i][j]
                
            rate = self.rowsum[i] / sumfrow
            residual = 0.0
            for j in range(self.ncol):
                fvalue = self.A[i][j] * rate
                ivalue = int(fvalue)
                if residual + fvalue - ivalue >= 0.5:
                    ivalue += 1
                residual += fvalue - ivalue
                image[i][j] = ivalue
                
        self.intA = imat
        return
            
    # after fIPF, integerize the seed values, while keeping the sum of them to match rowsum and colsum
    def iIPF(self):
        main_iteration = self.fIPF(self.A, self.rowsum, self.colsum)

        self.intA, self.fA = self.real_to_integer_and_float(self.A)
        # keep the row and col sum before modify the min and max values
        frowsum, fcolsum = self.marginals(self.fA)

        prev_sumfA = 0.0
        while 1:
            self.drop_minvalues()
            sub_iteration = self.fIPF(self.fA, frowsum, fcolsum)

            updated = self.rollup_maxvalues(frowsum, fcolsum)
            sumfA, sumfrow, sumfcol = self.nonezerocells_to_ones(frowsum, fcolsum)
            if abs(sumfA) < 1.0e-8 or (sumfrow == 0 and sumfcol==0) : break
            if abs(prev_sumfA - sumfA) > 0.01:
                if self.report:
                    print sumfA, sumfrow, sumfcol
            prev_sumfA = sumfA


if __name__ == '__main__':

    ## seed = [[204, 0, 90, 13, 0], [214, 0, 94, 13, 0], [160, 0, 48, 1, 0], [125, 0, 42, 3, 0],
    ##         [181, 0, 163, 0, 0], [71, 0, 119, 0, 0], [453, 0, 109, 0, 0], [193, 0, 30, 0, 0],
    ##         [171, 0, 49, 7, 0],  [447, 13, 241, 0, 0],  [500, 0, 76, 3, 0],  [155, 0, 7, 0, 0],
    ##         [148, 23, 211, 36, 0], [127, 51, 209, 39, 0], [6, 0, 8, 2, 1], [611, 3, 34, 0, 0],
    ##         [138, 5, 253, 63, 47], [516, 15, 82, 17, 0], [221, 0, 127, 0, 0], [0, 0, 1, 0, 0],
    ##         [45, 0, 26, 0, 0], [1, 0, 0, 0, 0]]
    ## rowsum = [316, 332, 264, 202, 373, 229, 576, 181, 258, 758, 601, 222, 471, 380, 18, 781, 483, 597, 512, 1, 103, 1]
    ## colsum = [4678, 138, 2535, 248, 60]

    seed = [[92, 44, 38, 8, 0],  [342, 0, 0, 0, 0], [438, 15, 28, 0, 0], [166, 0, 0, 0, 0],
            [394, 0, 0, 0, 0], [376, 19, 0, 0, 0], [402, 0, 0, 0, 0], [229, 0, 40, 0, 0],
            [482, 105, 29, 0, 0], [967, 273, 0, 0, 0], [352, 134, 243, 0, 0], [2, 1, 2, 0, 0],
            [353, 32, 232, 0, 0], [111, 490, 554, 0, 0], [231, 141, 13, 0, 0]]
    rowsum = [186,  428, 453, 183, 424, 427, 314, 271, 635, 1236, 731, 0, 647, 1034, 467]
    colsum = [4915, 1295, 1218, 8, 0]

    ipf = IPF_2dInteger(seed, rowsum, colsum)
    if ipf is None:
        print 'something is wrong'
    else:
        ipf.iIPF()
        print ipf.intA

    print 'done'
