#include <iostream>
#include <fstream>

using namespace std;

ifstream ifs;
ofstream ofs, ofsbw;

char r,g,b,h;

int main()
{
    ifs.open("lena.bmp");
    ofs.open("lena2.bmp");
    ofsbw.open("lena3.bmp");
	
	for(int i=0; i<54; i++)
	{
			ifs.get(h);
			ofs<<h;
			ofsbw<<h;
	}
    
    while(!ifs.eof())
    {
        ifs.get(r);
        ofs<<r;
        ifs.get(g);
        ofs<<g;
        ifs.get(b);
        ofs<<b;
       
        int lum = (int)( ((int)r + (int)g + (int)b)/3 );
			
        ofsbw<<(char)lum;
        ofsbw<<(char)lum;
        ofsbw<<(char)lum;
    }

    ifs.close();
    ofs.close();
    ofsbw.close();
    
    return 0;
}
