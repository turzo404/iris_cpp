#include<bits/stdc++.h>
using namespace std;

double sl[150], sw[150], pl[150], pw[150];

static void irisLoad()
{
    ifstream ip;
    ofstream op;

    string id, sL, sW, pL, pW, sp;

    ip.open("iris2.csv");
    op.open("iris3.txt");
    if(!ip.is_open())
        std::cout << "ERROR: File Open" << '\n';

    int i = 0;
    while(i<150)
    {
        getline(ip,id,',');
        getline(ip,sL,',');
        getline(ip,sW,',');
        getline(ip,pL,',');
        getline(ip,pW,',');
        getline(ip,sp,'\n');

        //std::cout << i;
        //std::cout <<    id  <<',';
        //std::cout <<    sL  <<',';
        //std::cout <<    sW  <<',';
        //std::cout <<    pL  <<',';
        //std::cout <<    pW  <<',';
        //std::cout <<    sp  <<'\n';

        sl[i]   =   atof(sL.c_str());
        sw[i]   =   atof(sW.c_str());
        pl[i]   =   atof(pL.c_str());
        pw[i]   =   atof(pW.c_str());

        op <<   sl[i]   <<',';
        op<<    sw[i]   <<',';
        op<<    pl[i]   <<',';
        op<<    pw[i]   <<',';
        op<<    sp      <<'\n';

        i++;
    }

    ip.close();
    op.close();
}

static void printArray()
{
    for(int i = 0; i<150; i++)
        cout<<sl[i]<<','<<sw[i]<<','<<pl[i]<<','<<pw[i]<<endl;
}

double sentosa_mean_sepalLength, sentosa_mean_sepalWidth;
double sentosa_mean_petalLength, sentosa_mean_petalWidth;
double virginica_mean_sepalLength, virginica_mean_sepalWidth;
double virginica_mean_petalLength, virginica_mean_petalWidth;
double versi_mean_sepalLength, versi_mean_sepalWidth;
double versi_mean_petalLength, versi_mean_petalWidth;

double sentosa_sd_sepalLength, sentosa_sd_sepalWidth;
double sentosa_sd_petalLength, sentosa_sd_petalWidth;
double virginica_sd_sepalLength, virginica_sd_sepalWidth;
double virginica_sd_petalLength, virginica_sd_petalWidth;
double versi_sd_sepalLength, versi_sd_sepalWidth;
double versi_sd_petalLength, versi_sd_petalWidth;

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
    cout<<"\n Sentosa : "<<endl;
    cout<<sentosa_mean_sepalLength<<" "<<sentosa_mean_sepalWidth<<" ";
    cout<<sentosa_mean_petalLength<<" "<<sentosa_mean_petalWidth<<" ";
    cout<<sentosa_sd_sepalLength<<" "<<sentosa_sd_sepalWidth<<" ";
    cout<<sentosa_sd_petalLength<<" "<<sentosa_sd_petalWidth<<" ";

    sentosa_sepalLength_left = sentosa_mean_sepalLength - sentosa_sd_sepalLength;
    sentosa_sepalLength_right = sentosa_mean_sepalLength + sentosa_sd_sepalLength;
    sentosa_sepalWidth_left =  sentosa_mean_sepalWidth - sentosa_sd_sepalWidth;
    sentosa_sepalWidth_right = sentosa_mean_sepalWidth + sentosa_sd_sepalWidth;
    sentosa_petalLength_left = sentosa_mean_petalLength - sentosa_sd_petalLength;
    sentosa_petalLength_right = sentosa_mean_petalLength + sentosa_sd_petalLength;
    sentosa_petalWidth_left = sentosa_mean_petalWidth - sentosa_sd_petalWidth;
    sentosa_petalWidth_right = sentosa_mean_petalWidth + sentosa_sd_petalWidth;


    cout<<"\n Versi : "<<endl;
    cout<<versi_mean_sepalLength<<" "<<versi_mean_sepalWidth<<" ";
    cout<<versi_mean_petalLength<<" "<<versi_mean_petalWidth<<" ";
    cout<<versi_sd_sepalLength<<" "<<versi_sd_sepalWidth<<" ";
    cout<<versi_sd_petalLength<<" "<<versi_sd_petalWidth<<" ";

    versi_sepalLength_left = versi_mean_sepalLength - versi_sd_sepalLength;
    versi_sepalLength_right = versi_mean_sepalLength + versi_sd_sepalLength;
    versi_sepalWidth_left =  versi_mean_sepalWidth - versi_sd_sepalWidth;
    versi_sepalWidth_right = versi_mean_sepalWidth + versi_sd_sepalWidth;
    versi_petalLength_left = versi_mean_petalLength - versi_sd_petalLength;
    versi_petalLength_right = versi_mean_petalLength + versi_sd_petalLength;
    versi_petalWidth_left = versi_mean_petalWidth - versi_sd_petalWidth;
    versi_petalWidth_right = versi_mean_petalWidth + versi_sd_petalWidth;

    cout<<"\n Virginica : "<<endl;
    cout<<virginica_mean_sepalLength<<" "<<virginica_mean_sepalWidth<<" ";
    cout<<virginica_mean_petalLength<<" "<<virginica_mean_petalWidth<<" ";
    cout<<virginica_sd_sepalLength<<" "<<virginica_sd_sepalWidth<<" ";
    cout<<virginica_sd_petalLength<<" "<<virginica_sd_petalWidth<<" ";

    virginica_sepalLength_left = virginica_mean_sepalLength - virginica_sd_sepalLength;
    virginica_sepalLength_right = virginica_mean_sepalLength + virginica_sd_sepalLength;
    virginica_sepalWidth_left =  virginica_mean_sepalWidth - virginica_sd_sepalWidth;
    virginica_sepalWidth_right = virginica_mean_sepalWidth + virginica_sd_sepalWidth;
    virginica_petalLength_left = virginica_mean_petalLength - virginica_sd_petalLength;
    virginica_petalLength_right = virginica_mean_petalLength + virginica_sd_petalLength;
    virginica_petalWidth_left = virginica_mean_petalWidth - virginica_sd_petalWidth;
    virginica_petalWidth_right = virginica_mean_petalWidth + virginica_sd_petalWidth;
}


int main()
{

    irisLoad();
    printArray();
    getMean();
    getSD();
    getRange();


    return 0;
}
