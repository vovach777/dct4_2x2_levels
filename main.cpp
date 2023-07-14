//vovach777 (c) 2023
//orinal work here:   https://godbolt.org/z/EabhKsqW1
//update:  https://godbolt.org/z/fGMKvn5qe
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <algorithm>



static const float _m = cosf(M_PI/4);


static void dct(float &a, float&b) {

    float aa = (a + b)/2.0f;
    float bb = a*_m - b*_m;
    a = aa; b = bb;
}

static void inv_dct(float &a, float&b) {

    float aa = a + b*_m;
    float bb = a - b*_m;
    a = aa; b = bb;
}

static void dct(int &a, int &b, int &c, int &d)
{
    float a_=a,b_=b,c_=c,d_=d;
	dct(a_,b_);
	dct(c_,d_);
	dct(a_,c_);
	dct(b_,d_);
    a = roundf(a_);
    b = roundf(b_);
    c = roundf(c_);
    d = roundf(d_);
}

static void inv_dct(int &a, int &b, int &c, int &d)
{
    float a_=a,b_=b,c_=c,d_=d;
	inv_dct(a_,c_);
	inv_dct(b_,d_);
	inv_dct(a_,b_);
	inv_dct(c_,d_);
    a = roundf(a_);
    b = roundf(b_);
    c = roundf(c_);
    d = roundf(d_);
}


using the_matrix = std::vector<std::vector<int>>;

std::ostream& operator << (std::ostream& o, const std::vector<std::vector<int>>& a)
{
    o << std::endl;
    for (auto & row : a)
    {
        for (auto v : row)
            o << std::setw(3) << v << " ";
        o << std::endl;
    }
    o << std::endl;
    return o;
}

the_matrix make_sample(int size) {
   the_matrix data(size,std::vector<int>(size));
  for (int y=0; y<size; y++)
  for (int x=0; x<size; x++)
    {
        double xx =  x * 2.1 - size;
        double yy =  y * 2.1 - size;
        double r = sqrt( xx*xx+yy*yy );
        //data[i] = r > 8.0 ? 0 : (int)r * 16 / 8;
        data[y][x] = (int)r*256/size/2.1;
        //data[i] = i % 8;
    }
    return std::move(data);
}


template <typename T>
static T median(T a, T b, T c) {
   if (a<b) {
      if (b<c)
         return b;
      return  (a < c) ? c : a;
   } else {
      if (b>=c)
         return b;
      return a < c ? a : c;
   }
}

void transform(the_matrix & m, int levels) {
    const int height=m.size();
    const int width=m[0].size();
    
    std::vector<std::array<int,3>> ring;
    std::array<int,3> top,left,topleft;
    
    //std::cout << width << "x" << height << std::endl;
    for (int level=1; level<=levels; ++level)
    {
        //std::cout << "level" << level << " " << ( 1<<(level) ) <<  " " << (1<<((level-1))) << std::endl;
                    
        const int Lwidth = width >> level;
        const int modR = Lwidth*2;
        ring.clear();
        ring.resize(modR);
                
                
        for (int h=0, r=modR, bsize=1<<(level), step=1<<((level-1)); h+step < height; h+=bsize)        
        {
            
        for (int x=0; x+step < width; x+=bsize,++r) {
            auto& b = m[h][x+step];
            auto& c = m[h+step][x];
            auto& d = m[h+step][x+step];
            dct( m[h][x], b, c, d);
            ring[r % modR] = std::array<int,3>{b,c,d};
            top   = h ? ring[(r - Lwidth) % modR ] : std::array<int,3>{0,0,0};
            left  = x ? ring[(r - 1 )  % modR ] : std::array<int,3>{0,0,0};
            topleft = x && h ? ring[(r - Lwidth - 1) % modR ] : std::array<int,3>{0,0,0};
            
            b -= median(top[0],left[0],top[0]+left[0]-topleft[0]);
            c -= median(top[1],left[1],top[1]+left[1]-topleft[1]);
            d -= median(top[2],left[2],top[2]+left[2]-topleft[2]);                    
            //std::cout << ",   " << (r-modR) << "=" << x;
        }
          //std::cout << std::endl << "r:"        << r  <<  " " << Lwidth << std::endl;
          //std::cout << std::endl;
        }
        //std::cout << "*** level " << level << " ***" << m;
    }
}


void inv_transform(the_matrix & m, int levels) {
    const int height=m.size();
    const int width=m[0].size();
    std::vector<std::array<int,3>> ring;
    std::array<int,3> top,left,topleft;
    //std::cout << width << "x" << height << std::endl;
    for (int level=levels; level>0; --level)
    {
        const int Lwidth = width >> level;
        const int modR = Lwidth*2;
        ring.clear();
        ring.resize(modR);
                
                
        for (int h=0, r=modR, bsize=1<<(level), step=1<<((level-1)); h+step < height; h+=bsize)        
        {
            
            for (int x=0; x+step < width; x+=bsize,++r) {
                top   = h ? ring[(r - Lwidth) % modR ] : std::array<int,3>{0,0,0};
                left  = x ? ring[(r - 1 )  % modR ] : std::array<int,3>{0,0,0};
                topleft = x && h ? ring[(r - Lwidth - 1) % modR ] : std::array<int,3>{0,0,0};
                
                auto& b = m[h][x+step];
                auto& c = m[h+step][x];
                auto& d = m[h+step][x+step];

                b += median(top[0],left[0],top[0]+left[0]-topleft[0]);
                c += median(top[1],left[1],top[1]+left[1]-topleft[1]);
                d += median(top[2],left[2],top[2]+left[2]-topleft[2]);                    

                ring[r % modR] = std::array<int,3>{b,c,d};

                inv_dct( m[h][x], b, c, d);
                
                //std::cout << ",   " << (r-modR) << "=" << x;
            }
          //std::cout << std::endl << "r:"        << r  <<  " " << Lwidth << std::endl;
          //std::cout << std::endl;
        }

    }
}


int main() {

  auto data = make_sample(16);
  std::cout << "matrix: " << data;
  transform(data,4);
  std::cout << "transformed: " << data;
  inv_transform(data,4);
  std::cout << "reconstructed: " << data;

}
