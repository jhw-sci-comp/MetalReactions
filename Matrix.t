#ifndef MATRIX_T
#define MATRIX_T

#include <vector>
#include <iostream>



template <typename T> class Matrix {
  public:
    // Constructors
    Matrix() = default;
	Matrix(unsigned int rows, unsigned int cols);
	Matrix(Matrix<T>&);
	// Destructor
	~Matrix();
	// Getter
	const unsigned int getRows() const;
	const unsigned int getCols() const;
	const std::vector<T>* getEntries() const;
	
		
	const T getEntry(unsigned int, unsigned int) const;
	void setEntry(unsigned int, unsigned int, T);
	const std::vector<T> getRow(unsigned int) const;
	void setRow(unsigned int, unsigned int, std::vector<T>);
	void swapRows(unsigned int, unsigned int);	
	void transpose();
	void createSubmatrix(unsigned int, unsigned int, Matrix<T>&) const;

	void setCol(int, std::vector<T>);
	void setCol(int, T);
	
	std::vector<T> operator* (const std::vector<T>&) const;
	Matrix<T>      operator* (const Matrix<T>&) const;
	Matrix<T>      operator* (const T) const;
	Matrix<T>      operator+ (const Matrix<T>&) const;
	Matrix<T>      operator- (const Matrix<T>&) const;
	Matrix<T>&     operator= (const Matrix<T>&);
		  
  private:
	unsigned int rows = 0;
	unsigned int cols = 0;
	std::vector<T>* entries = nullptr;
};




template <typename T>
Matrix<T> operator*(T scalar, const Matrix<T>& m) {
  Matrix<T> sm_product(m.getRows(), m.getCols());

  for (unsigned int i {0}; i < m.getRows(); i++) {    
    for (unsigned int j {0}; j < m.getCols(); j++) {
	  sm_product.setEntry(i, j, scalar * m.getEntry(i, j));
	}	
  }
  
  return sm_product;
}


template <typename T>
T det(const Matrix<T>& m) {
  T det_val = 0;
  T sign = 1;
  Matrix<T> submatrix(m.getRows() - 1, m.getCols() - 1);  

  
  if(m.getRows() == 1 && m.getCols() == 1) {	
	det_val = m.getEntry(0,0);
  }
  else if(m.getRows() == 2 && m.getCols() == 2) {	
	det_val = m.getEntry(0,0) * m.getEntry(1,1) - m.getEntry(1,0) * m.getEntry(0,1);
  }
  else {	
	for(unsigned int j = 0; j <= m.getCols() - 1; j++) {	  
	  m.createSubmatrix(0,j,submatrix);	  
	  det_val += sign * m.getEntry(0, j) * det(submatrix);	  
	  sign = -sign;
	}
  } 

  return det_val;
}



template <typename T>
Matrix<T>::Matrix( unsigned int rows_, unsigned int cols_) 
: rows {rows_}, cols {cols_} {	
  this->entries = new std::vector<T>(rows * cols, 0);
}

template <typename T>
Matrix<T>::Matrix(Matrix<T>& m) {
	this->rows = m.getRows();
	this->cols = m.getCols();

	this->entries = new std::vector<T>(rows * cols, 0);

	for(unsigned int i=0; i <= m.getRows() - 1; i++) {
		for(unsigned int j=0; j <= m.getCols() - 1; j++) {
			this->setEntry(i,j,getEntry(i,j));
		}
	}
}

template <typename T>
Matrix<T>::~Matrix() {	
	rows = 0;
	cols = 0;
  delete entries;	
}

template <typename T>
const unsigned int Matrix<T>::getRows() const {
  return rows;
}

template <typename T>
const unsigned int Matrix<T>::getCols() const {
  return cols;
}

template <typename T>
const std::vector<T>* Matrix<T>::getEntries() const {
  return entries;
}


template <typename T>
const T Matrix<T>::getEntry(unsigned int i, unsigned int j) const {  
  T entry = entries->at(i * cols + j);
  
  return entry;
}

template <typename T>
void Matrix<T>::setEntry(unsigned int rowIndex, unsigned int colIndex, T val) {
  this->entries->at(rowIndex * this->cols + colIndex) = val;
}


template <typename T>
const std::vector<T> Matrix<T>::getRow(unsigned int rowIndex) const {
  std::vector<T> row_temp;

  for(unsigned int j = 0; j < this->cols; j++) {
    row_temp.push_back(this->entries->at(rowIndex * this->cols + j));
  }

  return row_temp;
}


template <typename T>
void Matrix<T>::setRow(unsigned int rowIndex, unsigned int colIndex, std::vector<T> newRow) {
	for(unsigned int j = 0; j < colIndex; j++) {
		this->setEntry(rowIndex, j, newRow.at(j));
	}
}


template <typename T>
void Matrix<T>::swapRows(unsigned int rowIndex1, unsigned int rowIndex2) {
  T entry;

  for(unsigned int j = 0; j < this->cols; j++) {
    entry = this->entries->at(rowIndex1 * this->cols + j);
    this->entries->at(rowIndex1 * this->cols + j) = this->entries->at(rowIndex2 * this->cols + j);
    this->entries->at(rowIndex2 * this->cols + j) = entry;
  }
}




template <typename T>
void Matrix<T>::transpose() {    //read columns write rows
  unsigned int rows_temp {rows};
  std::vector<T>* temp_entries = new std::vector<T>(rows * cols); 
    
    
  for (unsigned int i {0}; i < entries->size(); i++) {
    temp_entries->at(i) = entries->at(i);
  }
  
  for(unsigned int i {0}; i < rows; i++) {
    for(unsigned int j {0}; j < cols; j++) {	  
	  entries->at(j * rows + i) = temp_entries->at(i * cols + j);
	}
  }
  
  rows = cols;
  cols = rows_temp; 
  
  delete temp_entries;  
}



template <typename T>
void Matrix<T>::createSubmatrix(unsigned int rowIndex, unsigned int colIndex, Matrix<T>& sub_matrix) const {
  Matrix<T> sub_temp(this->rows - 1, this->cols - 1);
  unsigned int i_sub {0}, j_sub {0};

  for(unsigned int i = 0; i <= this->rows - 2; i++) {
	if(i >= rowIndex) {
	  i_sub = i + 1;
	}
	else {
	  i_sub = i;
	}

	for(unsigned int j = 0; j <= this->cols - 2; j++) {
	  if(j >= colIndex) {
	    j_sub = j + 1;
	  }
	  else {
	    j_sub = j;
	  }			  

	  sub_temp.setEntry(i, j, this->getEntry(i_sub, j_sub));
	}
  }

  sub_matrix = sub_temp;  
}


template <typename T>
void Matrix<T>::setCol(int colIndex, std::vector<T> col_entries) {
	for(unsigned int i = 0; i < this->rows; i++) {
		this->setEntry(i, colIndex, col_entries.at(i));
	}
}


template <typename T>
void Matrix<T>::setCol(int colIndex, T val) {
	for(unsigned int i = 0; i < this->rows; i++) {
		this->setEntry(i, colIndex, val);
	}
}




template <typename T>
std::vector<T> Matrix<T>::operator* (const std::vector<T>& v) const {
  std::vector<T> mv_product;
  T sum;
  
  for (unsigned int i {0}; i < rows; i++) {
    sum = 0;
    for (unsigned int j {0}; j < cols; j++) {	
	  sum += entries->at(i * cols + j) * v.at(j) ;
	  
	}
	mv_product.push_back(sum);
  }
  
  return mv_product;
}


template <typename T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& m) const {
  Matrix<T> mm_product(rows, m.getCols());  
  T sum;
  
  
  for (unsigned int i {0}; i < rows; i++) {    
    for (unsigned int j {0}; j < m.getCols(); j++) {
	  sum = 0;
      for (unsigned int k {0}; k < cols; k++) {	  
		sum += this->getEntry(i, k) * m.getEntry(k, j);		
	  }	  
      mm_product.setEntry(i,j,sum);
	}	
  }
  
  return mm_product;
}


template <typename T>
Matrix<T>  Matrix<T>::operator* (const T scalar) const {  
  Matrix<T> sm_product(rows, cols);
  
  for (unsigned int i {0}; i < rows; i++) {    
    for (unsigned int j {0}; j < cols; j++) {
	  sm_product.setEntry(i,j, scalar * this->getEntry(i, j));
	}	
  }
  
  return sm_product;
}


template <typename T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& m) const {
  Matrix<T> mm_sum(rows, cols);  
  
  for (unsigned int i {0}; i < rows; i++) {    
    for (unsigned int j {0}; j < m.getCols(); j++) {
	  mm_sum.setEntry(i,j, this->getEntry(i, j) + m.getEntry(i, j));
	}	
  }
  
  return mm_sum;
}


template <typename T>
Matrix<T> Matrix<T>::operator- (const Matrix<T>& m) const {
  Matrix<T> mm_difference(rows, cols);  
  
  for (unsigned int i {0}; i < rows; i++) {    
    for (unsigned int j {0}; j < m.getCols(); j++) {
	  mm_difference.setEntry(i,j, this->getEntry(i, j) - m.getEntry(i, j));
	}	
  }
  
  return mm_difference;
}



template <typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& m) {
  this->rows = m.getRows();
  this->cols = m.getCols();
  
  for (unsigned int i {0}; i < rows; i++) {    
    for (unsigned int j {0}; j < m.getCols(); j++) {	  
	  this->setEntry(i, j, m.getEntry(i, j));
	}
  }

  return *this;
}


#endif