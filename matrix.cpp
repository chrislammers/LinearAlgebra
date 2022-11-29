#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class Matrix {
private:
  string name;
  int row;
  int col;
  int size;
  float * data;
public:

  //init
  Matrix() {
    name = "M";
    row = 1;
    col = 1;
    size = row*col;
    data = new float[1];
    data[0] = 0;
  }

  // init w/ some args
  Matrix(string n, int r, int c, float d[]) {
    name = n;
    row = r;
    col = c;
    size = row*col;
    data = d;
  }

  //init (only one number)
  Matrix(string n, int r, int c, float d) {
    name = n;
    row = r;
    col = c;
    size = row*col;
    float* array = new float[size];
    for(int jj = 0; jj < size; jj++) {
      // cout << "test " << jj << endl;
      array[jj] = d;
      // cout << array[jj] << endl;
    }
    data = array;
  }

  // //Elementary n*n matrix
  Matrix(int r){
    name = "I" + to_string(r);
    row = r;
    col = r;
    size = r*r;
    float* array = new float[r*r];
    int jj = 0;
    array[0] = 1;
    for(int ii = 1; ii < r*r-1; ++ii){
      // cout << ii << endl;
      for (int jj = 0; jj < r; ++jj){
        ++ii;
        array[ii] = 0;
        // cout << ii << endl;
      }
      array[ii] = 1;
    }
    data = array;
  }


  //Scalar multiplication
  void ScalarMult(float k) {
    for(int ii = 0; ii < this->size; ++ii) {
      this->data[ii] *= k;
    }
  }

  // Check if the arrays are the same size
  bool SameSize(Matrix m) {
    return (this->row == m.row && this->col == m.col);
  }

  Matrix operator+(const Matrix& m) {
    if(SameSize(m)) {
      string name = this->name + "+" + m.name;
      int r = m.row;
      int c = m.col;
      float* d = new float[r*c];
      for(int ii = 0; ii < r*c; ++ii) {
        d[ii] = this->data[ii] + m.data[ii];
      }
      return Matrix(name, r, c, d);
    }
    float zero = 0;
    cerr << this->name << "+" << m.name << " is an undefined matrix." << endl;
    return Matrix(this->name+m.name, 1, 1, zero);
  }

  // Matrix Multiplication!    Mat1 Columns = Mat2 Rows
  friend Matrix MatMult(Matrix a, Matrix b) {
    if(a.col == b.row) {
      string name = a.name + b.name;
      int r = a.row;
      int c = b.col;
      float* d = new float[r*c];
      float entry;
      for(int rr = 0; rr < r; ++rr) {
        for(int cc = 0 ; cc < c; ++cc) {
          entry = 0;
          // indexing a 2 dimensional matrix:
          //    location on line = current col + (totcal col x current row)
          for(int ii = 0; ii < a.col; ++ii) {
            // cout << a.data[ii+a.col*rr] << "\t" << b.data[ii+b.row*cc] <<endl;
            entry += a.data[ii+a.col*rr] * b.data[ii+b.row*cc];
          }
          // cout << endl;

          // cout << a.data[rr] << "*" << b.data[cc] << endl;
          // entry = cc+c*rr;
          // entry = a.data[rr] * b.data[cc];
          d[cc+c*rr] = entry;
        }
      }

      Matrix m = Matrix(name, r, c, d);
      return m;
    }
    cerr << a.name+b.name << " is an undefined matrix." << endl;
    float zero = 0;
    return Matrix(a.name+b.name, 1, 1, zero);
  }

  //Transpose: mirror matrix on the diagonal
  friend Matrix Transpose(Matrix a) {
    int r = a.col;
    int c = a.row;
    float* d = new float[a.size];
    for(int ii = 0; ii < a.col; ++ii) {
      for(int jj = 0; jj < a.row; ++jj) {
        // indexing a 2 dimensional matrix:
        //    location on line = current col + (total col x current row)
        d[jj+a.row*ii] = a.data[ii+a.col*jj];
      }
    }

    Matrix m = Matrix(a.name + "T", r, c, d);
    return m;
  }

  // Row Reduction (augmenting two matrices) with form [this | a]
  void RowReduce(Matrix a) {
    if(SameSize(a)) {
      a.RowDiv(0,index(0,0));
      RowDiv(0,index(0,0));
      // iterates through each row
      for(int ii = 0; ii < this->row; ++ii) {
        for(int jj = 0; jj < ii; ++jj) {
          cout << "ii = " << ii << "\tjj = " << jj << endl;
          a.RowSub(ii, jj, index(ii, jj), index(ii-1, jj));
          RowSub(ii, jj, index(ii, jj), index(ii-1, jj));
        }
        // a.RowDiv(ii, index(ii, ii));
        // RowDiv(ii, index(ii, ii));
      }
    } else {
      cerr << "I cannot augment matrices of two different sizes!" << endl;
    }
  }


  // row "x" = (a) row "x" -  (b) row "y"
  void RowSub(int x, int y, int a = 1, int b = 1){
    if(x != y) {
      // ii represents the location on each row
      for(int ii = 0; ii<this->col; ++ii) {
        this->data[ii + this->row*x] = a*(this->data[ii + this->row*x]) - b*(this->data[ii + this->row*y]);
      }
    } else {
      cerr << "I cannot subtract a row from itself!" << endl;
    }
  }

  // row "x" divided by float "n"
  void RowDiv(int x, float n){
    if (n != 0) {
      // ii represents the location on row "x"
      for(int ii = 0; ii<this->col; ++ii) {
        this->data[ii + this->row*x] /= n;
      }
    } else {
      cerr << "I cannot divide a row by 0" << endl;
    }
  }

  // Find number at position r, c
  float index(int row, int col) {
    return this->data[col + this->row*row];
  }

  // Print function
  void printM(int n) {
    // cout << this->data << endl;
    int count = 0;
    cout << this->name << " =\t" << this->row << " x " << this->col << endl;
    for(int ii = 0; ii < this->row; ++ii){
      cout << "\t[ ";
      for(int ii = 0; ii < this->col; ++ii) {
        // %.*f uses n as precision
        printf("%.*f", n, this->data[count]);
        cout << " ";
        count += 1;
      }
      cout << "]\n";
    }
  }
};

int main(){
  // cout << "Test" << endl;
  Matrix m1 = Matrix();
  float d2[] = {2,4,6,8,10,12};
  float d4[] = {1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41};
  Matrix m2 = Matrix("M2", 2, 3, d2);
  Matrix m3 = Matrix("M3", 3, 7, 4);
  Matrix m4 = Matrix("M4", 3, 7, d4);
  // m2.ScalarMult(4);
  m1.printM(2);
  m2.printM(2);
  m3.printM(2);
  m4.printM(2);

  Matrix m2m4 = MatMult(m2, m4);
  m2m4.printM(2);

  Matrix m5 = m3+m4;
  m5.printM(2);
  m5.ScalarMult(0.5);
  m5.printM(2);

  Matrix i7 = Matrix(7);
  i7.printM(1);

  Matrix m6 = Transpose(m5);
  m6.printM(2);

  Matrix m6m4 = MatMult(m6, m4);
  m6m4.ScalarMult(0.01);
  m6m4.printM(1);

  // m6m4.RowDiv(2, 5);
  // m6m4.RowSub(3,-1);
  m6m4.RowReduce(i7);
  m6m4.printM(1);
  i7.printM(1);


  // cout << m6m4.index(2,5);

  return 0;
}
