///////////////////////////////////////////////////////////////////////////////////
///
///
///		Naive Bayesian Classification on Iris Data Set
///		Turzo Ahsan Sami
///		Roll: 181812, MIT-18
///
///
///////////////////////////////////////////////////////////////////////////////////



#include<bits/stdc++.h>
using namespace std;


//////////////////////////////////////////////////////////////////////
//
//	reading iris data from file
//	and storing into separate arrays based upon the attributes
//
//////////////////////////////////////////////////////////////////////

#define dataSize 150

double sl[150], sw[150], pl[150], pw[150];

static void irisLoad()
{
    ifstream ip;

    string id, sL, sW, pL, pW, sp;

    ip.open("iris2.csv");
    if(!ip.is_open())
        std::cout << "ERROR: File Open" << '\n';

    int i = 0;
    while(!ip.eof())
    {
        getline(ip,id,',');
        getline(ip,sL,',');
        getline(ip,sW,',');
        getline(ip,pL,',');
        getline(ip,pW,',');
        getline(ip,sp,'\n');

        sl[i]   =   atof(sL.c_str());
        sw[i]   =   atof(sW.c_str());
        pl[i]   =   atof(pL.c_str());
        pw[i]   =   atof(pW.c_str());

        i++;
    }

    ip.close();
}

static void printArray()
{
     for(int i = 0; i<dataSize; i++)
         cout<<sl[i]<<" , "<<sw[i]<<" , "<<pl[i]<<" , "<<pw[i]<<endl;
}



//////////////////////////////////////////////////////////////////////
//
//	calculating Mean and Standard Deviation of each attribute
//
//////////////////////////////////////////////////////////////////////




double sentosa_mean_sepalLength, sentosa_mean_sepalWidth;
double sentosa_mean_petalLength, sentosa_mean_petalWidth;
double sentosa_sd_sepalLength, sentosa_sd_sepalWidth;
double sentosa_sd_petalLength, sentosa_sd_petalWidth;

double versi_mean_sepalLength, versi_mean_sepalWidth;
double versi_mean_petalLength, versi_mean_petalWidth;
double versi_sd_sepalLength, versi_sd_sepalWidth;
double versi_sd_petalLength, versi_sd_petalWidth;

double virginica_mean_sepalLength, virginica_mean_sepalWidth;
double virginica_mean_petalLength, virginica_mean_petalWidth;
double virginica_sd_sepalLength, virginica_sd_sepalWidth;
double virginica_sd_petalLength, virginica_sd_petalWidth;


static void getMean()
{
    // Sentosa
    double t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    for(int i=0; i<50; i++)
    {
        t1 += sl[i];
        t2 += sw[i];
        t3 += pl[i];
        t4 += pw[i];
    }
    sentosa_mean_sepalLength    = t1/50;
    sentosa_mean_sepalWidth     = t2/50;
    sentosa_mean_petalLength    = t3/50;
    sentosa_mean_petalWidth     = t4/50;

    // Versi
    t1=t2=t3=t4=0;
    for(int i=50; i<100; i++)
    {
        t1 += sl[i];
        t2 += sw[i];
        t3 += pl[i];
        t4 += pw[i];
    }
    versi_mean_sepalLength    = t1/50;
    versi_mean_sepalWidth     = t2/50;
    versi_mean_petalLength    = t3/50;
    versi_mean_petalWidth     = t4/50;

    // Virginica
    t1=t2=t3=t4=0;
    for(int i=100; i<150; i++)
    {
        t1 += sl[i];
        t2 += sw[i];
        t3 += pl[i];
        t4 += pw[i];
    }
    virginica_mean_sepalLength    = t1/50;
    virginica_mean_sepalWidth     = t2/50;
    virginica_mean_petalLength    = t3/50;
    virginica_mean_petalWidth     = t4/50;

}

static void getSD()
{
    double t1, t2, t3, t4;

    // sentosa
    t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    for(int i=0; i<50; i++)
    {
        t1 += pow((sl[i] - sentosa_mean_sepalLength), 2);
        t2 += pow((sw[i] - sentosa_mean_sepalWidth), 2);
        t3 += pow((pl[i] - sentosa_mean_petalLength), 2);
        t4 += pow((pw[i] - sentosa_mean_petalWidth), 2);
    }
    sentosa_sd_sepalLength  = sqrt(t1/50);
    sentosa_sd_sepalWidth   = sqrt(t2/50);
    sentosa_sd_petalLength  = sqrt(t3/50);
    sentosa_sd_petalWidth   = sqrt(t4/50);

    // versi
    t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    for(int i=50; i<100; i++)
    {
        t1 += pow((sl[i] - versi_mean_sepalLength), 2);
        t2 += pow((sw[i] - versi_mean_sepalWidth), 2);
        t3 += pow((pl[i] - versi_mean_petalLength), 2);
        t4 += pow((pw[i] - versi_mean_petalWidth), 2);
    }
    versi_sd_sepalLength    = sqrt(t1/50);
    versi_sd_sepalWidth     = sqrt(t2/50);
    versi_sd_petalLength    = sqrt(t3/50);
    versi_sd_petalWidth     = sqrt(t4/50);

    // virginica
    t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    for(int i=100; i<150; i++)
    {
        t1 += pow((sl[i] - virginica_mean_sepalLength), 2);
        t2 += pow((sw[i] - virginica_mean_sepalWidth), 2);
        t3 += pow((pl[i] - virginica_mean_petalLength), 2);
        t4 += pow((pw[i] - virginica_mean_petalWidth), 2);
    }
    virginica_sd_sepalLength    = sqrt(t1/50);
    virginica_sd_sepalWidth     = sqrt(t2/50);
    virginica_sd_petalLength    = sqrt(t3/50);
    virginica_sd_petalWidth     = sqrt(t4/50);

}


//////////////////////////////////////////////////////////////////////////////////
//
//	calculating the range of each attribute : x = { (mean - sd), (mean + sd) }
//
/////////////////////////////////////////////////////////////////////////////////


double sentosa_sepalLength_left, sentosa_sepalLength_right;
double sentosa_sepalWidth_left, sentosa_sepalWidth_right;
double sentosa_petalLength_left, sentosa_petalLength_right;
double sentosa_petalWidth_left, sentosa_petalWidth_right;

double versi_sepalLength_left, versi_sepalLength_right;
double versi_sepalWidth_left, versi_sepalWidth_right;
double versi_petalLength_left, versi_petalLength_right;
double versi_petalWidth_left, versi_petalWidth_right;

double virginica_sepalLength_left, virginica_sepalLength_right;
double virginica_sepalWidth_left, virginica_sepalWidth_right;
double virginica_petalLength_left, virginica_petalLength_right;
double virginica_petalWidth_left, virginica_petalWidth_right;

static void getRange()
{
    sentosa_sepalLength_left = sentosa_mean_sepalLength - sentosa_sd_sepalLength;
    sentosa_sepalLength_right = sentosa_mean_sepalLength + sentosa_sd_sepalLength;
    sentosa_sepalWidth_left =  sentosa_mean_sepalWidth - sentosa_sd_sepalWidth;
    sentosa_sepalWidth_right = sentosa_mean_sepalWidth + sentosa_sd_sepalWidth;
    sentosa_petalLength_left = sentosa_mean_petalLength - sentosa_sd_petalLength;
    sentosa_petalLength_right = sentosa_mean_petalLength + sentosa_sd_petalLength;
    sentosa_petalWidth_left = sentosa_mean_petalWidth - sentosa_sd_petalWidth;
    sentosa_petalWidth_right = sentosa_mean_petalWidth + sentosa_sd_petalWidth;

    versi_sepalLength_left = versi_mean_sepalLength - versi_sd_sepalLength;
    versi_sepalLength_right = versi_mean_sepalLength + versi_sd_sepalLength;
    versi_sepalWidth_left =  versi_mean_sepalWidth - versi_sd_sepalWidth;
    versi_sepalWidth_right = versi_mean_sepalWidth + versi_sd_sepalWidth;
    versi_petalLength_left = versi_mean_petalLength - versi_sd_petalLength;
    versi_petalLength_right = versi_mean_petalLength + versi_sd_petalLength;
    versi_petalWidth_left = versi_mean_petalWidth - versi_sd_petalWidth;
    versi_petalWidth_right = versi_mean_petalWidth + versi_sd_petalWidth;

    virginica_sepalLength_left = virginica_mean_sepalLength - virginica_sd_sepalLength;
    virginica_sepalLength_right = virginica_mean_sepalLength + virginica_sd_sepalLength;
    virginica_sepalWidth_left =  virginica_mean_sepalWidth - virginica_sd_sepalWidth;
    virginica_sepalWidth_right = virginica_mean_sepalWidth + virginica_sd_sepalWidth;
    virginica_petalLength_left = virginica_mean_petalLength - virginica_sd_petalLength;
    virginica_petalLength_right = virginica_mean_petalLength + virginica_sd_petalLength;
    virginica_petalWidth_left = virginica_mean_petalWidth - virginica_sd_petalWidth;
    virginica_petalWidth_right = virginica_mean_petalWidth + virginica_sd_petalWidth;
}

int sentosaSepalLengthCount = 0;
int sentosaSepalWidthCount = 0;
int sentosaPetalLengthCount = 0;
int sentosaPetalWidthCount = 0;

int versiSepalLengthCount = 0;
int versiSepalWidthCount = 0;
int versiPetalLengthCount = 0;
int versiPetalWidthCount = 0;

int virginicaSepalLengthCount = 0;
int virginicaSepalWidthCount = 0;
int virginicaPetalLengthCount = 0;
int virginicaPetalWidthCount = 0;

double sepalLength[3];
double sepalWidth[3]; 
double petalLength[3];
double petalWidth[3];

static void createRange()
{
	for(int i=0; i<3; i++)
		sepalLength[i]=sepalWidth[i]=petalLength[i]=petalWidth[i]=0;
	
	double max=-999, min=999;
	for(int i=0; i<dataSize; i++)
	{
		if(max<=sl[i])	max=sl[i];
		if(min<=sl[i])	min=sl[i];
	}
	sepalLength[0] = min - max;
	sepalLength[1] = (min+max)/2;
	sepalLength[2] = max + min;
	
	max=-999, min=999;
	for(int i=0; i<dataSize; i++)
	{
		if(max<=sW[i])	max=sW[i];
		if(min<=sW[i])	min=sW[i];
	}
	sepalWidth[0] = 0;
	sepalWidth[1] = (min + max)/2;
	sepalWidth[2] = max + max;
	
	max=-999, min=999;
	for(int i=0; i<dataSize; i++)
	{
		if(max<=pl[i])	max=pl[i];
		if(min<=pl[i])	min=pl[i];
	}
	petalLength[0] = min - max;
	petalLength[1] = (min+max)/2;
	petalLength[2] = max + min;
	
	max=-999, min=999;
	for(int i=0; i<dataSize; i++)
	{
		if(max<=pw[i])	max=pw[i];
		if(min<=pw[i])	min=pw[i];
	}
	petalWidth[0] = min - max;
	petalWidth[1] = (min+max)/2;
	petalWidth[2] = max + min;
}

static void count()
{
	for(int i = 0; i<dataSize; i++)
	{
		if(sl[i] >= sentosa_sepalLength_left && sl[i]<= sentosa_sepalLength_right)
			sentosaSepalLengthCount++;
		if(sl[i] >= versi_sepalLength_left && sl[i]<= versi_sepalLength_right)
			versiSepalLengthCount++;
		if(sl[i] >= virginica_sepalLength_left && sl[i]<= virginica_sepalLength_right)
			virginicaSepalLengthCount++;
	}
	
	for(int i = 0; i<dataSize; i++)
	{
		if(sw[i] >= sentosa_sepalWidth_left && sw[i]<= sentosa_sepalWidth_right)
			sentosaSepalWidthCount++;
		if(sw[i] >= versi_sepalWidth_left && sw[i]<= versi_sepalWidth_right)
			versiSepalWidthCount++;
		if(sw[i] >= virginica_sepalWidth_left && sw[i]<= virginica_sepalWidth_right)
			virginicaSepalWidthCount++;
	}
	
	for(int i = 0; i<dataSize; i++)
	{
		if(pl[i] >= sentosa_petalLength_left && pl[i]<= sentosa_petalLength_right)
			sentosaPetalLengthCount++;
		if(pl[i] >= versi_petalLength_left && pl[i]<= versi_petalLength_right)
			versiPetalLengthCount++;
		if(pl[i] >= virginica_petalLength_left && pl[i]<= virginica_petalLength_right)
			virginicaPetalLengthCount++;
	}
	
	for(int i = 0; i<dataSize; i++)
	{
		if(sl[i] >= sentosa_petalWidth_left && sl[i]<= sentosa_petalWidth_right)
			sentosaPetalWidthCount++;
		if(sl[i] >= versi_petalWidth_left && sl[i]<= versi_petalWidth_right)
			versiPetalWidthCount++;
		if(sl[i] >= virginica_petalWidth_left && sl[i]<= virginica_petalWidth_right)
			virginicaPetalWidthCount++;
	}
}

int  SL[3]; SL[0] = SL[1] = SL[2] = 0;
int  SW[3]; SW[0] = SW[1] = SW[2] = 0;
int  PL[3]; PL[0] = PL[1] = PL[2] = 0;
int  PW[3]; PW[0] = PW[1] = PW[2] = 0;

static void populate()
{
	for(int i=0; i<dataSize; i++)
	{
		if(sl[i]>=sepalLength[0] && sl[i]<=sepalLength[0])
			SL[0]++;
		if(sl[i]>=sepalLength[1] && sl[i]<=sepalLength[1])
			SL[1]++;
		if(sl[i]>=sepalLength[2] && sl[i]<=sepalLength[2])
			SL[2]++;
			
		if(sw[i]>=sepalWidth[0] && sw[i]<=sepalWidth[0])
			SW[0]++;
		if(sw[i]>=sepalWidth[1] && sw[i]<=sepalWidth[1])
			SW[1]++;
		if(sw[i]>=sepalWidth[2] && sw[i]<=sepalWidth[2])
			SW[2]++;
			
		if(pl[i]>=petalLength[0] && pl[i]<=petalLength[0])
			PL[0]++;
		if(pl[i]>=petalLength[1] && pl[i]<=petalLength[1])
			PL[1]++;
		if(pl[i]>=petalLength[2] && pl[i]<=petalLength[2])
			PL[2]++;
			
		if(pw[i]>=petalWidth[0] && pw[i]<=petalWidth[0])
			PW[0]++;
		if(pw[i]>=petalWidth[1] && pw[i]<=petalWidth[1])
			PW[1]++;
		if(pw[i]>=petalWidth[2] && pw[i]<=petalWidth[2])
			PW[2]++;
	}
}

double cenP = 0;

static void isSentosa(double sll, double sww, double pll, double pww)
{
	if(sll>=SL[0] && sll<SL[1])
	{
		if(sww>=SW[0] && sww<SW[1])
		{
			if(pll>=PL[0] && pll<PL[1])
			{
				if(pww>=PW[0] && pww<PW[1])
			}
		}
	}
	else if(sll>=SL[1] && sll<SL[2])
	{
		if(pll>=)
	}						
	else
	{
		
	}	
}



///////////////////////////////////////////////////////////////////////////////////
///
///
///		main Function
///
///
///////////////////////////////////////////////////////////////////////////////////



int main()
{

    irisLoad();
    printArray();
    getMean();
    getSD();
    getRange();
	createRange();
	count();
	populate();
	
    double sll, sww, pll, pww;
    cin>>sll>>sww>>pll>>pww;
	
	isSentosa(sll, sww, pll, pww);
	
    return 0;
}
