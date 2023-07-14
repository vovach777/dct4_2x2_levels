//vovach777 (c) 2023
//orinal work here:   https://godbolt.org/z/EabhKsqW1
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <algorithm>
#ifdef INT_DCT4
void dct4(int &a, int&b, int&c, int&d) {

    int aa = (a + b + c + d)/4;
    int bb = a + b*2/5 - c*2/5 - d;
    int cc = a*3/4 - b*3/4 - c*3/4 + d*3/4;
    int dd = a*2/5 - b + c - d*2/5;
    a = aa; b = bb; c = cc; d = dd;
}

void idct4(int &a, int&b, int&c, int&d) {

    int aa = a*2 + b + c*3/4 + d*2/5;
    int bb = a*2 + b*2/5 - c*3/4 - d;
    int cc = a*2 - b*2/5 - c*3/4 + d;
    int dd = a*2 - b + c*3/4 - d*2/5; 
    a = aa/2; b = bb/2; c = cc/2; d = dd/2;
}
#endif

/*
0.92388 0.382683 -0.382683 -0.92388 
0.707107 -0.707107 -0.707107 0.707107 
0.382683 -0.92388 0.92388 -0.382683 
*/

    constexpr float _h = 0.92388f;
    constexpr float _m = 0.707107f;
    constexpr float _l = 0.382683f;

static void dct4f(int &a, int&b, int&c, int&d) {

    float aa = (a + b + c + d)/4.0f;
    float bb = a*_h + b*_l - c*_l - d*_h;
    float cc = a*_m - b*_m - c*_m + d*_m;
    float dd = a*_l - b*_h + c*_h - d*_l;
    a = round(aa); b = round(bb); c = round(cc); d = round(dd);
}

static void idct4f(int &a, int&b, int&c, int&d) {

    float aa = a*2.0f + b*_h + c*_m + d*_l;
    float bb = a*2.0f + b*_l - c*_m - d*_h;
    float cc = a*2.0f - b*_l - c*_m + d*_h;
    float dd = a*2.0f - b*_h + c*_m - d*_l; 
    a = round(aa/2.0f); b = round(bb/2.0f); c = round(cc/2.0f); d = round(dd/2.0f);
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

static int get_m_value(const the_matrix &m, int x, int y)
{
    if ( x < 0 || x >= m[0].size() || y < 0 || y >= m.size())
        return 0;
    else
        return m[y][x];
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
            dct4f( m[h][x], b, c, d);
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

                idct4f( m[h][x], b, c, d);
                
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
