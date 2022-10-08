class kll {
private:
  vector<vector<double>> kllArray;
  int height;
  int n;
public:
  kll() {}

  void init(int _n, vector<double> w) {
    kllArray.push_back(w);
    n = _n;
    int lastLevel = 0;
    while (_n > 2) {
      vector<double> tmpArray;
      int lastSize = kllArray[lastLevel].size();
      for (int i = 0;i < lastSize;i += 2) {
        double tmp = kllArray[lastLevel][i];
        if (lastSize > i + 1) tmp += kllArray[lastLevel][i + 1];
        tmpArray.push_back(tmp);
      }
      kllArray.push_back(tmpArray);
      lastLevel++;
      _n = ceil(double(_n) / 2);
    }
    height = lastLevel;
  }

  ~kll() {}

  int query(double preSum) {
    int tmpLevel = height;
    //the highest Level can contain one value or two values. we have to check one by one.
    int highestLen = kllArray[tmpLevel].size();
    int i;
    for (i = 0;i < highestLen;i++) {
      if (preSum > kllArray[tmpLevel][i]) preSum -= kllArray[tmpLevel][i];
      else break;
    }
    while (tmpLevel > 0) {
      int leftSon = i * 2, rightSon = i * 2 + 1;
      tmpLevel--;
      if (kllArray[tmpLevel][leftSon] > preSum) {
        i = leftSon;
      }
      else {
        i = rightSon;
        preSum -= kllArray[tmpLevel][leftSon];
      }
    }
    return i;
  }

  /**
   * kllArray[0][index] = 0 by default. If not, you should delete it first. kllArray[0][index]=value after  `insert` operator.
   */
  void modifyValue(double value, int idx) {
    kllArray[0][idx] = value;
    int tmpLevel = 1;
    while (tmpLevel <= height) {
      idx /= 2;
      kllArray[tmpLevel][idx] += value;
      tmpLevel++;
    }
  }

  void addOperator(double value) {
    n++;
    int _n = n;
    kllArray[0].push_back(value);
    int idx = n - 1;
    int tmpLevel = 1;
    while (tmpLevel <= height) {
      idx /= 2;
      if (kllArray[tmpLevel].size() <= idx) kllArray[tmpLevel].push_back(value);
      else kllArray[tmpLevel][idx] += value;
      tmpLevel++;
      _n = ceil(double(_n) / 2);
    }
    if (_n > 2) {
      vector<double> tmpArray;
      int lastSize = kllArray[height].size();
      for (int i = 0;i < lastSize;i += 2) {
        double tmp = kllArray[height][i];
        if (lastSize > i + 1) tmp += kllArray[height][i + 1];
        tmpArray.push_back(tmp);
      }
      kllArray.push_back(tmpArray);
      height++;
    }
  }

  double resetValue(int idx) {
    double value = kllArray[0][idx];
    kllArray[0][idx] = 0;
    if (idx == n - 1) kllArray[0].pop_back();
    int tmpLevel = 1;
    while (tmpLevel <= height) {
      idx /= 2;
      kllArray[tmpLevel][idx] -= value;
      if (kllArray[tmpLevel][idx] == 0 && kllArray[tmpLevel].size() == idx + 1) kllArray[tmpLevel].pop_back();
      tmpLevel++;
    }
    if (kllArray[height].size() < 2) {
      kllArray.pop_back();
      height--;
    }
    return value;
  }

  void deleteOperator(int idx) {
    if (idx >= n) {
      cout << "delete out of the bound" << endl;
      exit(-1);
    }
    resetValue(idx);
    if (idx == n - 1) {
      n--;
      return;
    }
    double value = resetValue(n - 1);
    modifyValue(value, idx);
    n--;
  }
};