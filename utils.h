#include <string.h>
unsigned long long peak_mem()
{
    char buf[1024];
    FILE *fp;
    unsigned long long peaksize;
    fp = fopen("/proc/self/status", "r");
    if (fp == NULL)
    {
        cout << "fp = null" << endl;
        return -1;
    }
    while (fgets(buf, sizeof(buf) - 1, fp) != NULL)
    {
        if (sscanf(buf, "VmPeak:%llu", &peaksize) > 0)
        {
            break;
        }
    }
    fclose(fp);
    return peaksize;
}

int mkpath(string s, mode_t mode = 0755)
{
    size_t pre = 0, pos;
    string dir;
    int mdret;
    if (s[s.size() - 1] != '/')
    {
        s += '/';
    }
    while ((pos = s.find_first_of('/', pre)) != string::npos)
    {
        dir = s.substr(0, pos++);
        pre = pos;
        if (dir.size() == 0)
            continue;
        if ((mdret = ::mkdir(dir.c_str(), mode)) && errno != EEXIST)
        {
            return mdret;
        }
    }
    return mdret;
}

void readGraphAttr(string filelabel, string filedir, uint &n, uint &m)
{
    string graphAttr = filedir + filelabel + ".attribute";
    cout << "Read graph attributes..." << endl;
    cout << "FilePath: " << graphAttr.c_str() << endl;
    ifstream graphAttrIn(graphAttr.c_str());
    string tmp;
    graphAttrIn >> tmp >> n;
    graphAttrIn >> tmp >> m;
    cout << "n=" << n << endl;
    graphAttrIn.close();
    return;
}

void readGraph(string filelabel, string filedir, uint *outEL, uint *outPL, double *outWEL, uint n, uint m, int sorted)
{

    string neiNode, neiWeight, neiNum;

    if (sorted)
    {
        neiNode = filedir + filelabel + ".outSortedEdges";
        neiWeight = filedir + filelabel + ".outSortedWEdges";
        neiNum = filedir + filelabel + ".outPtr";
    }
    else
    {
        neiNode = filedir + filelabel + ".outEdges";
        neiWeight = filedir + filelabel + ".outWEdges";
        neiNum = filedir + filelabel + ".outPtr";
    }

    cout << "Read graph ..." << endl;

    ifstream neiNumIn(neiNum.c_str(), ios::in | ios::binary);
    ifstream neiWeightIn(neiWeight.c_str(), ios::in | ios::binary);
    ifstream neiNodeIn(neiNode.c_str(), ios::in | ios::binary);
    neiNodeIn.read((char *)&outEL[0], sizeof(outEL[0]) * m);
    neiWeightIn.read((char *)&outWEL[0], sizeof(outWEL[0]) * m);
    neiNumIn.read((char *)&outPL[0], sizeof(outPL[0]) * (n + 1));
    neiNumIn.close();
    neiWeightIn.close();
    neiNodeIn.close();
    return;
}

void argParser(int argc, char **argv, string &filedir, string &filelabel, long &querynum, vector<double> &epss, uint &L)
// void argParser(int argc, char **argv, string &filedir, string &filelabel, long &querynum, double &eps, uint &L)
{
    char *endptr;
    int i = 1;
    while (i < argc)
    {
        if (!strcmp(argv[i], "-d"))
        {
            if (++i < argc)
                filedir = argv[i];
        }
        else if (!strcmp(argv[i], "-f"))
        {
            if (++i < argc)
                filelabel = argv[i];
        }
        else if (!strcmp(argv[i], "-e"))
        {
            while (++i < argc && argv[i][0] != '-')
            {
                // eps = strtod(argv[++i], &endptr);
                double eps = strtod(argv[i], &endptr);
                if ((eps == 0 || eps > 1) && endptr)
                {
                    cerr << "Invalid eps argument" << endl;
                    exit(1);
                }
                epss.push_back(eps);
            }
            i--;
        }
        else if (!strcmp(argv[i], "-qn"))
        {
            if (++i < argc)
                querynum = strtod(argv[i], &endptr);
            if ((querynum < -2) && endptr)
            {
                cerr << "Invalid querynum argument:" << querynum << endl;
                exit(1);
            }
        }
        else if (!strcmp(argv[i], "-l"))
        {
            int tempL;
            if (++i < argc)
                tempL = int(strtod(argv[i], &endptr));
            if ((tempL < 0) && endptr)
            {
                cerr << "Invalid hop number" << endl;
                exit(1);
            }
            L = uint(tempL);
        }
        else
        {
            cerr << "Err args:" << argv[i] << endl;
            exit(1);
        }
        i++;
    }
}