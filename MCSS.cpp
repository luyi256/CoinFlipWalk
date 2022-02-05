class Graph
{
public:
    uint n; //number of nodes
    uint m; //number of edges
    string filedir;
    string filelabel;
    double *totw;
    uint *out_size_set;
    unordered_map<uint, node *> *umap;       //<neighbor_id, nodeptr> umap[source_id]
    unordered_map<uint, node *> *chain_head; //<low_bound,chain_head_ptr> chain_head[source_id]
    Graph(string _filedir, string _filelabel);
    void add(uint s, uint t, double w);
    void modify(uint s, uint t, double w);
    void del(uint s, uint t);
    ~Graph()
    {
    }
};
struct node
{
    node *pre;
    node *post;
    uint idx;
    double w;
};

void Graph::add(uint s, uint t, double w)
{
    unordered_map<uint, node *>::const_iterator got;
    got = umap[s].find(t);
    if (got != umap[s].end())
        return; //已存在
    node *post_node, *cur_node;
    int listnum = ceil(log2(w));
    if (listnum <= 0)
        listnum = 0;
    got = chain_head[s].find(listnum);
    if (got != chain_head[s].end())
    {
        post_node = got->second;
        cur_node = new node{
            NULL,
            post_node,
            t,
            w};
        post_node->pre = cur_node;
    }
    else
    {
        cur_node = new node{
            NULL,
            NULL,
            t,
            w};
    }
    chain_head[s][listnum] = cur_node;
    umap[s][t] = cur_node;
    out_size_set[s]++;
    totw[s] += w;
}

void Graph::modify(uint s, uint t, double w)
{
    unordered_map<uint, node *>::const_iterator got;
    got = umap[s].find(t);
    if (got == umap[s].end())
        return; //不存在
    totw[s] += w - umap[s][t]->w;
    umap[s][t]->w = w;
}

void Graph::del(uint s, uint t)
{
    unordered_map<uint, node *>::const_iterator got;
    got = umap[s].find(t);
    if (got == umap[s].end())
        return; //不存在
    node *cur_node;
    cur_node = got->second;
    if (cur_node->pre == NULL)
    {
        int listnum = ceil(log2(cur_node->w));
        if (listnum <= 0)
            listnum = 0;
        node *head = cur_node->post;
        chain_head[s][listnum] = head;
        if (head != NULL)
            head->pre = NULL;
    }
    else
    {
        cur_node->pre->post = cur_node->post;
        if (cur_node->post != NULL)
            cur_node->post->pre = cur_node->pre;
    }

    totw[s] -= cur_node->w;
    out_size_set[s]--;
    free(cur_node);
}

Graph::Graph(string _filedir, string _filelabel)
{
    filedir = _filedir;
    filelabel = _filelabel;
    stringstream ss_dir, ss_attr, ss_outEL, ss_outPL, ss_outWEL, ss_inEL, ss_inPL, ss_inWEL;
    uint *outEL;
    uint *outPL;
    double *outWEL;
    ss_attr << filedir << filelabel << ".attribute";
    ifstream in_attr;
    in_attr.open(ss_attr.str());
    cout << "FilePath: " << ss_attr.str() << endl;
    cout << "Read graph attributes..." << endl;
    string tmp;
    in_attr >> tmp >> n;
    in_attr >> tmp >> m;
    cout << "n=" << n << " m=" << m << endl;

    in_attr.close();

    cout << "Read graph edges..." << endl;

    ss_outEL << filedir << filelabel << ".outEdges";
    ss_outWEL << filedir << filelabel << ".outWEdges";

    ss_outPL << filedir << filelabel << ".outPtr";

    outEL = new uint[m];
    outPL = new uint[n + 1];
    outWEL = new double[m];
    //每个节点的出节点下标从哪里开始
    ifstream outpf(ss_outPL.str(), ios::in | ios::binary);
    outpf.read((char *)&outPL[0], sizeof(outPL[0]) * (n + 1));
    //所有节点的出节点编号
    ifstream outf, outef;
    outf.open(ss_outEL.str(), ios::in | ios::binary);
    outf.read((char *)&outEL[0], sizeof(outEL[0]) * m);
    outef.open(ss_outWEL.str(), ios::in | ios::binary);
    outef.read((char *)&outWEL[0], sizeof(outWEL[0]) * m);
    outf.close();
    outpf.close();
    outef.close();
    chain_head = new unordered_map<uint, node *>[n];
    umap = new unordered_map<uint, node *>[n];
    totw = new double[n];
    out_size_set = new uint[n];
    for (int i = 0; i < n; i++)
    {
        node *post_node, *cur_node;
        uint outSize = (outPL[i + 1] - outPL[i]);
        out_size_set[i] = outSize;
        totw[i] = 0;
        unordered_map<uint, node *>::const_iterator got;
        for (int j = 0; j < outSize; j++)
        {
            uint idx = outEL[(outPL[i] + j)];
            double w = outWEL[(outPL[i] + j)];
            totw[i] += w;
            int listnum = ceil(log2(w));
            if (listnum <= 0)
                listnum = 0;
            got = chain_head[i].find(listnum);
            if (got != chain_head[i].end())
            {
                post_node = got->second;
                cur_node = new node{
                    NULL,
                    post_node,
                    idx,
                    w};
                post_node->pre = cur_node;
            }
            else
            {
                cur_node = new node{
                    NULL,
                    NULL,
                    idx,
                    w};
            }
            chain_head[i][listnum] = cur_node;
            umap[i][idx] = cur_node;
        }
    }
}