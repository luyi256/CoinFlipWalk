double max(double a, double b)
{
    return a > b ? a : b;
}
struct AVLnode
{
    double value;
    AVLnode* left, * right;
    int leftHeight, rightHeight;
    int leftCount, rightCount;
    double leftSum, rightSum;
    AVLnode(double v)
    {
        value = v;
        left = NULL;
        right = NULL;
        leftHeight = 0;
        rightHeight = 0;
        leftCount = 1;
        rightCount = 0;
        leftSum = value;
        rightSum = 0;
    }

    void setValue(double v)
    {
        value = v;
        if (left)
            leftSum = left->leftSum + left->rightSum + value;
        else
            leftSum = value;
    }

    void setLeft(AVLnode* node)
    {
        left = node;
        if (node)
        {
            leftHeight = max(node->leftHeight, node->rightHeight) + 1;
            leftCount = node->leftCount + node->rightCount + 1;
            leftSum = node->leftSum + node->rightSum + value;
        }
        else
        {
            leftHeight = 0;
            leftCount = 1;
            leftSum = value;
        }
    }

    void setRight(AVLnode* node)
    {
        right = node;
        if (node)
        {
            rightHeight = max(node->leftHeight, node->rightHeight) + 1;
            rightCount = node->leftCount + node->rightCount;
            rightSum = node->leftSum + node->rightSum;
        }
        else
        {
            rightHeight = 0;
            rightCount = 0;
            rightSum = 0;
        }
    }

    AVLnode* rotateLeft()
    {
        AVLnode* b = right;
        setRight(b->left);
        b->setLeft(this);
        return b;
    }

    AVLnode* rotateRight()
    {
        AVLnode* b = left;
        setLeft(b->right);
        b->setRight(this);
        return b;
    }

    int factor()
    {
        return rightHeight - leftHeight;
    }
};

AVLnode* addNode(AVLnode* root, uint index, AVLnode* node)
{
    if (root == NULL)
        return node;
    if (index < root->leftCount)
    {
        root->setLeft(addNode(root->left, index, node));
        if (root->factor() < -1)
        {
            if (root->left->factor() > 0)
                root->setLeft(root->left->rotateLeft());
            return root->rotateRight();
        }
    }
    else
    {
        root->setRight(addNode(root->right, index - root->leftCount, node));
        if (root->factor() > 1)
        {
            if (root->right->factor() < 0)
                root->setRight(root->right->rotateRight());
            return root->rotateLeft();
        }
    }
    return root;
}

AVLnode* updateNode(AVLnode* root, uint index, double v)
{
    if (root == NULL)
        return root;
    if (index + 1 < root->leftCount)
    {
        root->setLeft(updateNode(root->left, index, v));
    }
    else if (index + 1 > root->leftCount)
    {
        root->setRight(updateNode(root->right, index - root->leftCount, v));
    }
    else
        root->setValue(v);
    return root;
}

uint prefixSumIndex(AVLnode* root, double sum, uint indexAdd)
{
    if (root == NULL)
        return 0;
    if (root->left == NULL && root->right == NULL) {//leaf
        if (root->value >= sum) return root->leftCount + indexAdd;
        else return root->leftCount + indexAdd + 1;
    }
    if (root->leftSum < sum)
    {
        return prefixSumIndex(root->right, sum - root->leftSum, indexAdd + root->leftCount);
    }
    else // sum<=leftSum
    {
        if (root->left == NULL) //just the root
            return indexAdd + 1;
        else
            return prefixSumIndex(root->left, sum, indexAdd);
    }
}

AVLnode* deleteLast(AVLnode* root)
{
    if (root == NULL)
        return NULL;
    if (root->right)
    {
        root->setRight(deleteLast(root->right));
        if (root->factor() < -1)
        {
            if (root->left->factor() > 0)
                root->setLeft(root->left->rotateLeft());
            return root->rotateRight();
        }
        return root;
    }
    else if (root->left)
    {
        AVLnode* tmp = root->left;
        delete root;
        return tmp;
    }
    else
    {
        delete root;
        return NULL;
    }
}