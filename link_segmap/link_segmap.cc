#include <set>
#include <map>
#include <list>
#include <vector>
#include <stdexcept>
#include <complex>
#include <fitsio.h>
#include <iostream>

typedef unsigned long ulong;
typedef unsigned int uint;

// helper class: Grid for a image
// provides mapping between pixel numbers and coordinates
class Grid {
 public:
  Grid() {}
  Grid(uint N0_, uint N1_) : N0(N0_), N1(N1_) {}
  uint operator() (ulong pixel, bool direction) const {
    if (direction)
      return pixel/N0;
    else
      return pixel%N0;
  }
  ulong getPixel(uint x, uint y) const {
    if (x >= 0 && x < N0 && y >= 0 && y < N1)
      return x + y*N0;
    else
      return -1;
  }
  ulong size() const {
    return N0*N1;
  }
  ulong getNeighborPixel(ulong pixel, unsigned char direction) const {
    uint x = operator()(pixel, 0), y = operator()(pixel, 1);
    ulong index;
    switch(direction) {
    case 0: 
      // the pixel itself
      index = pixel; // y*N0 + x;
      break;
    case 1: 
      if (x>0) index = y*N0 + x - 1; // left
      else index = -1;
      break;
    case 2:
      if (x<N0-1) index = y*N0 + x + 1;  // right neighbour
      else index = -1;
      break;
    case 3: 
      if (y>0 && x>0) index = (y-1)*N0 + x - 1;  // bottom left
      else index = -1;
      break;   
    case 4: 
      if (y>0) index = (y-1)*N0 + x;  // bottom
      else index = -1;
      break;
    case 5: 
      if (y>0 && x<N0-1) index = (y-1)*N0 + x + 1;  // bottom right
      else index = -1;
      break;  
    case 6: 
      if (y<N1-1 && x>0) index = (y+1)*N0 + x - 1;  // top left
      else index = -1;
      break;  
    case 7: 
      if (y<N1-1) index = (y+1)*N0 + x ;  // top
      else index = -1;
      break;
    case 8:
      if (y<N1-1 && x<N0-1) index = (y+1)*N0 + x + 1;  // top right
      else index = -1;
      break;
    }
    return index;
  }
 private:
  uint N0, N1;
};


// templated version of cfitsio datatypes
template <typename T> inline
int getDataType(const T& entry) {
  // default type, uses template specialization for other types
  // see below
  return TBYTE;
}
template<> inline
int getDataType<int>(const int& entry) {
  return TINT;
}
template<> inline
int getDataType<unsigned int>(const unsigned int& entry) {
  return TUINT;
}
template<> inline
int getDataType<long>(const long& entry) {
  return TLONG;
}
template<> inline
int getDataType<unsigned long>(const unsigned long& entry) {
  return TULONG;
}
template<> inline
int getDataType<float>(const float& entry) {
  return TFLOAT;
}
template<> inline
int getDataType<double>(const double& entry) {
  return TDOUBLE;
}
template<> inline
int getDataType<std::complex<float> >(const std::complex<float>& entry) {
  return TCOMPLEX;
}
template<> inline
int getDataType<std::complex<double> >(const std::complex<double>& entry) {
  return TDBLCOMPLEX;
}
template<> inline
int getDataType<std::string>(const std::string& entry) {
  return TSTRING;
}


// helper class
// stores the contents of FITS file in a std::vector
template <class T>
class Image : public std::vector<T> {
public:
  Image(std::string filename) {
    int status = 0, naxis;
    fitsfile* fptr;
    fits_open_file(&fptr, filename.c_str(), false, &status);
    if (status != 0)
      throw std::runtime_error("Image: Cannot open " + filename);
    fits_get_img_dim(fptr, &naxis, &status);
    if (naxis!=2)
      throw std::invalid_argument("Image: naxis != 2." + filename + " does not hold image");
    
    long naxes[2] = {1,1};
    fits_get_img_size(fptr, naxis, naxes, &status);
    grid = Grid(naxes[0],naxes[1]);
    std::vector<T>::resize(grid.size());
    long firstpix[2] = {1,1};
    T val;
    int datatype = getDataType(val);
    fits_read_pix(fptr, datatype, firstpix, grid.size(), NULL, &std::vector<T>::front(), NULL, &status);
    if (status != 0)
      throw std::runtime_error("Image: Cannot read FITS image from " + filename);
    fits_close_file(fptr, &status);
    if (status != 0)
      throw std::runtime_error("FITS: Cannot close FITS file " + filename);
  }
  Grid grid;
};


// Find a list of connected pixels which have a particular value
// This is a Friend-of-Friend algorithm with a linking length of 1 pixel.
// It starts by putting startpixel into the pixelset and adds all other
// connected pixels with the same value (ignoring those it has seen already).
// For performance reason, the segmap will be altered (flipping the sign
// of an objects segmap pixel if it has been included in a group).
template <class T>
std::list<unsigned long> linkPixels(Image<T>& segmap, T value, ulong startpixel) {
  std::list<ulong> pixellist;
  pixellist.push_back(startpixel);

  std::list<ulong>::iterator iter = pixellist.begin();
  uint pixelnumber = 0;
  const Grid& grid = segmap.grid;
  while (pixelnumber < pixellist.size()) {
    ulong pixel = *iter;
    // loop over all direct neighbors and check if they have same value
    for (uint dir = 1; dir <= 8 ; dir++) {
      long neighbor = grid.getNeighborPixel(pixel,dir);
      if (neighbor != -1) {
	if (segmap[neighbor] == value) {
	  pixellist.push_back(neighbor);
	  segmap[neighbor] *= -1; // remove them from unseen list
	}
      }
    }
    iter++;
    pixelnumber++;
  }
  return pixellist;
}


int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " <segmap> " << std::endl;
    exit(1);
  }
  
  // open segmap
  // need to allow negative numbers, therefore cast to signed long
  Image<long> segmap(argv[1]);
  std::multimap<long, std::list<ulong> > counter;
  for (ulong i = 0; i < segmap.size(); i++) {
    long id = segmap[i];
    if (id > 0) {
      std::list<ulong> pixellist = linkPixels(segmap, id, i);
      counter.insert(std::pair<long, std::list<ulong> >(id, pixellist));
    }
  }
  
  // simple output the ID and # of fragments
  std::cout << "# ID\tPIXELS\tX_CENTER\tY_CENTER" << std::endl;
  for (std::multimap<long, std::list<ulong> >::iterator iter = counter.begin(); iter != counter.end(); iter++) {
    std::cout << iter->first << "\t" << iter->second.size();
    double x = 0, y = 0;
    for (std::list<ulong>::iterator piter = iter->second.begin(); piter != iter->second.end(); piter++) {
      x += segmap.grid(*piter, 0);
      y += segmap.grid(*piter, 1);
    }
    x /= iter->second.size();
    y /= iter->second.size();
    std::cout << "\t" << x << "\t" << y << std::endl;
  }
}
