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
        if (p == 0)
            return 0;
        if (k > 1e-6)
            return p;
        return p - 1;
    }
    void add(double k) //在尾部加一个k的值
    {
        n++;
        BITree.resize(n + 1);
        BITree[n] = k;
        if (n % 2 == 1)
            return;
        uint mul = 1;
        while (n % (mul * 2) == 0)
        {
            mul *= 2;
            BITree[n] += BITree[n - mul / 2];
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
