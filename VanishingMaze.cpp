// C++11
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <queue>

using namespace std;

int n, startPos;
array<int, 4> ADJ;
vector<int> cells;
vector<bool> broken;
int current;
int pending;
int score;

class Maze
{
};

vector<int> bfs(int p0, int target)
{
    queue<int> q;
    q.push(p0);
    vector<int> prev((n + 2) * (n + 2), -2);
    prev[p0] = -1;
    vector<int> res;
    while (!q.empty())
    {
        int v = q.front();
        q.pop();
        // cerr << "target : " << target << " " << v << endl;
        for (int a : ADJ)
        {
            int t = v + a;
            if (cells[t] == 0 || broken[t] || prev[t] != -2)
            {
                continue;
            }

            prev[t] = v;
            if (cells[t] == target || cells[t] == -1)
            {
                for (int u = t; u != p0; u = prev[u])
                {
                    // cerr << "find : " << target << " " << u << endl;
                    res.push_back(u);
                }
                reverse(res.begin(), res.end());
                return res;
            }

            q.push(t);
        }
    }
    return res;
}

void simulate(vector<int> &pos)
{
    for (int p : pos)
    {
        if (cells[p] == current || cells[p] == -1)
        {
            cells[p] = current;
            current++;
            score += pending;
            pending = 0;
        }
        if (cells[p] > current)
        {
            broken[p] = true;
            pending += current * cells[p];
        }
    }
}

string toResult(vector<int> &moves)
{
    string s;
    cerr << moves.size() << endl;
    for (int i = 1; i < (int)moves.size(); i++)
    {
        int d = moves[i] - moves[i - 1];
        auto j = distance(ADJ.begin(), find(ADJ.begin(), ADJ.end(), d));
        s += "LURD"[j];
        cerr << moves[i] / (n + 2) << " " << moves[i] % (n + 2) << endl;
        cerr << d << " " << s << endl;
    }
    return s;
}

string greedy()
{
    int p = startPos;
    vector<int> resultMoves;
    resultMoves.push_back(p);
    for (int target = 1; target <= n; target++)
    {
        auto ps = bfs(p, target);
        if (ps.empty())
        {
            break;
        }
        for (int i = 0; i < (int)ps.size(); i++)
        {
            cerr << current << " " << ps[i] / (n + 2) << " " << ps[i] % (n + 2) << endl;
        }
        copy(ps.begin(), ps.end(), back_inserter(resultMoves));
        simulate(ps);
        p = ps.back();
    }
    return toResult(resultMoves);
}

class VanishingMaze
{
public:
    string getPath(vector<int> numbers, int playerRow, int playerCol, double power)
    {
        n = (int)sqrt(numbers.size());
        ADJ[0] = -1;
        ADJ[2] = 1;
        ADJ[1] = -(n + 2);
        ADJ[3] = (n + 2);
        cells.resize((n + 2) * (n + 2));
        fill(cells.begin(), cells.end(), 0);
        for (int y = 0; y < n; y++)
        {
            for (int x = 0; x < n; x++)
            {
                cells[(y + 1) * (n + 2) + (x + 1)] = numbers[y * n + x];
            }
        }
        startPos = (playerRow + 1) * (n + 2) + playerCol + 1;
        broken.assign((n + 2) * (n + 2), false);
        current = 1;

        auto res = greedy();
        cerr << "calcScore : " << score * pow((current - 1) * 1.0 / n, power) << endl;

        return res;
    }
};
// -------8<------- end of solution submitted to the website -------8<-------

template <class T>
void getVector(vector<T> &v)
{
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main()
{
    VanishingMaze vm;
    int S2;
    cin >> S2;
    vector<int> numbers(S2);
    getVector(numbers);
    int R, C;
    double P;
    cin >> R >> C >> P;

    string ret = vm.getPath(numbers, R, C, P);
    cout << ret << endl;
    cout.flush();
}
