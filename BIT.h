class BIT
{
public:
    vector<double> BITree;
    int n;
    int find(double k) //返回的是在数组中的位置，非第几个
    {
        int p = 0;
        int q = 1 << int(floor(log2(n)));
        while (q != 0)
        {
            if (p + q <= n)
            {
                double m = BITree[p + q];
                if (k > m || ((fabs)(k - m) < 1e-6))
                {
                    p = p + q;
                    k = k - m;
                }
            }
            q = q >> 1;
        }
        if (k > 0)
            return p;
        else
            return p - 1;
    }
    void resize(int _n)
    {
        if (_n > n)
        {
            BITree.resize(_n + 1);
            for (uint j = n + 1; j <= _n; j++)
                BITree[j] = 0;
            n = _n;
        }
    }

    void updateBIT(int index, double val) //在数组中位置是index，排行index+1
    {
        index = index + 1;
        while (index <= n)
        {
            BITree[index] += val;
            index += index & (-index);
        }
    }

    BIT(int _n, double arr[])
    {
        n = _n;
        BITree = vector<double>(n + 1, 0);
        for (int i = 0; i < n; i++)
            updateBIT(i, arr[i]);
    }
    BIT() {}
    ~BIT() {}
};
