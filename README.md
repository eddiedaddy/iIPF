# iIPF
A IPF (Iterative Proportional Fitting) class that matches controls of row or column, or both. 

Iterative Proportional Fitting (IPF) is a method commonly used in modeling, where vectors are expended to matrices.
For example in transportation modeling, trip origin in column vector and trip destination in row vector are used to create
an O/D matrix to be fed into traffic assignment step.  Few methods of converting the two vectors to a matrix are used, including
Fratar, and Gravity model.  Those are virtually same procedure by iteratively adjusting the cell values in the matrix to match
column sums and row sums.  Difference between Fratar and Gravity models is how the initial cell values are developed.  For Fratar,
1 and 0 are used depend on the geographic connectivity, while gravity model starts with the "impedance" of travel between the origin
and destination.

Socio-economic data development also heavily relies on IPF too.  In many cases, pieces of information from Census are combined to
create more meaningfull joint tables, like number of households by size by income, which is not the table available directly from
the Census, yet highly desired data in urban modeling.

The method is alternatively matching the sums of rows and columns as adjusting the cell values in the matrix.  Overall the matching
is improved as the iterations go on.  But in any case, it is hard, if possible, to match the all of row sums and column sum together,
especially if the controls (row margins and column margins) are in integers, as well as the expected result is also a matrix of
integer values

This class is to match the integer controls when it is done.  The first "i" represents "integer". The implemented idea is as following

- input :
  rowsum  (sum of rowsum = sum of colsum = total to be filled in the result matrix in integer)
  colsum
  seed (matrix)

- steps :
1) perform normal IPF with real numbers
2) separate the "real" matrix to "integer" and "fractional" matrices.
3) make new colsum and rowsum vectors from the "fractional" matrix.  Since the sums should be integers too.
4) disturbe the "fractional" matrix, by removing the smallest none-zero value
5) perform normal IPF with new row and column sums and disturbed fractional matrix.
6) check the largest number in the IPF'd fractional matrix.  
    If the largests one is more than 1,
    - add 1 to the integer matrix to the same cell location of the largest one from the fractional matrix
    - remove 1 from rowsum and colsum
    - set to 0 to the cell value of the largest one in the fractional matrix
7) check the remainder, if the sum of rowsum, colsum are zeros, stop.  Otherwise go back to step 4)

