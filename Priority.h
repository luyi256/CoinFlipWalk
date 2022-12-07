
class Priority
{
public:
  vector<int> heap;
  Priority(){}
  Priority(vector<int> sorted)
  {
    heap.resize(sorted.size());
    memcpy(&heap[0], &sorted[0], sizeof(int) * sorted.size());
  }

  ~Priority()
  {
  }

  void initPriority(vector<int> sorted)
  {
    heap.resize(sorted.size());
    memcpy(&heap[0], &sorted[0], sizeof(int) * sorted.size());
  }

  void addKey(double newKey)
  {
    heap.push_back(newKey);
    int father = heap.size() / 2;
    while (father > 0)
    {
      int leftsonidx = father * 2;
      int rightsonidx = father * 2 + 1;
      int maxson = leftsonidx;
      if (rightsonidx <= heap.size() && heap[rightsonidx - 1] > heap[leftsonidx - 1])
        maxson = rightsonidx;
      if (heap[maxson - 1] > heap[father - 1])
      {
        int tmp = heap[maxson - 1];
        heap[maxson - 1] = heap[father - 1];
        heap[father - 1] = tmp;
      }
      else
        break;
      father /= 2;
    }
  }
};