#undef _GLIBCXX_DEBUG

#include <algorithm>
#include <functional>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <numeric>
#include <cassert>
#include <sys/time.h>

using namespace std;

#ifdef LOCAL
const double TL = 2.5;
#else
const double TL = 10.0;
#endif

double get_time() {
  timeval tv; 
  gettimeofday(&tv, 0); 
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

struct rand_gen {
  static const int MAX = 2147483647;
  static constexpr double Q_MAX = 1.0 / MAX;

  int x = 8753, y = 239017, z = 1000000123;

  inline int next_int() {
    int t = x ^ (x << 11);
    x = y;
    y = z;
    z ^= (z >> 19) ^ t ^ (t >> 8);
    return z;
  }

  double next_double() {
    return next_int() * Q_MAX;
  }
} rng;

template<typename T>
void shuffle(vector<T> &a) {
  int n = (int) a.size();
  for (int i = 1; i < n; i++) {
    swap(a[i], a[rng.next_int() % (i + 1)]);
  }
}

string to_string(string s) {
  return '"' + s + '"';
}

string to_string(const char* s) {
  return to_string((string) s);
}

string to_string(char c) {
  string res = "";
  res += c;
  return res;
}

string to_string(bool b) {
  return (b ? "true" : "false");
}

template <typename A, typename B>
string to_string(pair<A, B> p) {
  return "(" + to_string(p.first) + ", " + to_string(p.second) + ")";
}

template <typename A>
string to_string(A v) {
  bool first = true;
  string res = "{";
  for (const auto &x : v) {
    if (!first) {
      res += ", ";
    }
    first = false;
    res += to_string(x);
  }
  res += "}";
  return res;
}

void debug_out() { cerr << endl; }

template <typename Head, typename... Tail>
void debug_out(Head H, Tail... T) {
  cerr << " " << to_string(H);
  debug_out(T...);
}

#ifdef LOCAL
#define debug(...) cerr << "[" << #__VA_ARGS__ << "]:", debug_out(__VA_ARGS__)
#else
#define debug(...) 42
#endif

class dsu {
  public:
  vector<int> p;
  int n;

  dsu(int _n) : n(_n) {
    p.resize(n);
    iota(p.begin(), p.end(), 0);
  }

  inline int get(int x) {
    return (x == p[x] ? x : (p[x] = get(p[x])));
  }

  inline bool unite(int x, int y) {
    x = get(x);
    y = get(y);
    if (x != y) {
      p[x] = y;
      return true;
    }
    return false;
  }
};

const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = {0, -1, 0, 1};
const string dz = "ULDR";

const int inf = (int) 1e9;

double start_time;

vector<int> numbers;
int start;
double power;
int n;
int s;
vector<vector<int>> nei;
vector<vector<int>> by_number; // by_number[i] : cell positions of the number i

string build_by_rm(const vector<int> &rm) {
  vector<int> board = numbers;
  int pos = start;
  int num = 1;
  string ret = "";
  for (int x : rm) {
    vector<int> que(1, pos);
    vector<int> was(n * n, 0);
    vector<int> pr(n * n, -1);
    was[pos] = 1;
    for (int b = 0; b < (int) que.size(); b++) {
      if (b > 0 && (board[que[b]] == 0 || board[que[b]] == -1 || board[que[b]] >= num)) {
        continue;
      }
      for (int dir = 0; dir < 4; dir++) {
        int to = nei[que[b]][dir];
        if (!was[to]) {
          was[to] = 1;
          pr[to] = dir;
          que.push_back(to);
        }
      }
    }
    assert(was[x] == 1);
    string ops = "";
    int id = x;
    while (id != pos) {
      ops += dz[pr[id]];
      id = nei[id][pr[id] ^ 2];
    }
    ret += string(ops.rbegin(), ops.rend());
    if (board[x] == -1 || board[x] == num) {
      board[x] = num;
      num++;
    } else {
      board[x] = 0;
    }
    pos = x;
  }
  return ret;
}

double calc(const vector<int> &rm) {
  if (!rm.empty() && rm[0] == start) {
    return -1;
  }
  vector<int> board = numbers;
  int num = 1;
  dsu d(n * n);
  int pos = start;
  int pending = 0;
  int ans = 0;
  for (int x : rm) {
    assert(num <= s);
    if (board[x] == 0) {
      return -1;
    }
    if (board[x] > 0 && board[x] < num) {
      return -1;
    }
    int ok = 0;
    for (int z : nei[pos]) {
      if (z == x) {
        ok = 1;
        break;
      }
      if (board[z] > 0 && board[z] < num) {
        for (int y : nei[x]) {
          if (d.get(y) == d.get(z)) {
            ok = 1;
            break;
          }
        }
      }
    }
    if (!ok) {
      return -1;
    }
    if (board[x] > num) {
      pending += board[x];
      board[x] = 0;
    } else {
      ans += num * pending;
      pending = 0;
      if (board[x] == -1) {
        board[x] = num;
        for (int y : nei[x]) {
          if (board[y] > 0 && board[y] <= num) {
            d.unite(y, x);
          }
        }
      }
      for (int z : by_number[num]) {
        if (board[z] == 0) {
          continue;
        }
        for (int y : nei[z]) {
          if (board[y] > 0 && board[y] <= num) {
            d.unite(y, z);
          }
        }
      }
      num++;
    }
    pos = x;
  }
  double score = ans;
  score *= pow((double) (num - 1) / s, power);
  return score;
}

vector<pair<int,int>> opt;

void optimize(vector<int> &rm, double &score) {
  vector<int> nums(rm.size());
  nums[0] = 1;
  for (int i = 0; i < (int) rm.size() - 1; i++) {
    nums[i + 1] = nums[i] + (numbers[rm[i]] == nums[i] || numbers[rm[i]] == -1);
  }
  vector<int> useful(n * n, 0);
  for (int x : rm) {
    useful[x] = 1;
  }
  vector<int> tries(n * n, 0);
  opt.clear();
  for (int i = (int) rm.size() - 1; i >= 0; i--) {
    vector<int> board = numbers;
    int nnum = 1;
    for (int ii = 0; ii < i; ii++) {
      int& bb = board[rm[ii]];
      if (bb == nnum || bb == -1) {
        bb = nnum++;
      } else {
        assert(bb > nnum);
        bb = 0;
      }
    }
    dsu dd(n * n);
    for (int j = 0; j < n * n; j++) {
      if (board[j] > 0 && board[j] < nums[i]) {
        for (int x : nei[j]) {
          if (board[x] > 0 && board[x] < nums[i]) {
            dd.unite(j, x);
          }
        }
      }
    }
    for (int j = 0; j < n * n; j++) {
      if (!useful[j] && numbers[j] > nums[i]) {
        if (tries[j] > 0) {
          continue;
        }
        int ok = 0;
        for (int x : nei[j]) {
          if (x == rm[i]) {
            ok = 1;
            break;
          }
          if (board[x] > 0 && board[x] < nums[i]) {
            for (int y : nei[rm[i]]) {
              if (dd.get(x) == dd.get(y)) {
                ok = 1;
                break;
              }
            }
          }
        }
        if (!ok) {
          continue;
        }
        int pr = (i == 0 ? start : rm[i - 1]);
        ok = 0;
        for (int x : nei[j]) {
          if (x == pr) {
            ok = 1;
            break;
          }
          if (board[x] > 0 && board[x] < nums[i]) {
            for (int y : nei[pr]) {
              if (dd.get(x) == dd.get(y)) {
                ok = 1;
                break;
              }
            }
          }
        }
        if (!ok) {
          continue;
        }
        tries[j]++;
        rm.insert(rm.begin() + i, j);
        double new_score = calc(rm);
        if (new_score > -0.5) {
          assert(new_score > score);
          score = new_score;
          useful[j] = 1;
          opt.emplace_back(i, j);
        } else {
          rm.erase(rm.begin() + i);
        }
      }
    }
  }
  string ret = build_by_rm(rm);
  if ((int) ret.size() > 4 * n * s) {
    return;
  }
  // optimize-new:
  vector<int> rm_before_new = rm;
  vector<pair<int,int>> opt_before_new = opt;
  for (int i = (int) rm.size() - 2; i >= 0; i--) {
    vector<int> que(1, rm[i]);
    vector<int> was(n * n, 0);
    vector<int> pr(n * n, -1);
    was[rm[i]] = 1;
    for (int b = 0; b < (int) que.size(); b++) {
      if (b > 0 && (useful[que[b]] || numbers[que[b]] <= nums.back())) {
        continue;
      }
      for (int dir = 0; dir < 4; dir++) {
        int to = nei[que[b]][dir];
        if (b == 0 && to == rm[i + 1]) {
          continue;
        }
        if (!was[to]) {
          was[to] = 1;
          pr[to] = dir;
          que.push_back(to);
        }
      }
    }
    if (!was[rm[i + 1]]) {
      continue;
    }
    int id = rm[i + 1];
    int added = 0;
    while (true) {
      id = nei[id][pr[id] ^ 2];
      if (id == rm[i]) {
        break;
      }
      useful[id] = 1;
      rm.insert(rm.begin() + i + 1, id);
      opt.emplace_back(i + 1, id);
      added++;
    }
//    debug(i, added, rm);
    double new_score = 1.0;
    if (new_score > -0.5) {
      i += added + 1;
    } else {
      assert(false);
      rm.erase(rm.begin() + i + 1, rm.begin() + i + added + 1);
    }
  }
  ret = build_by_rm(rm);
  if ((int) ret.size() > 4 * n * s) {
    rm = rm_before_new;
    opt = opt_before_new;
  }
  score = calc(rm);
}

class VanishingMaze {
public:
  string getPath(vector<int> _numbers, int playerRow, int playerCol, double _power) {
    start_time = get_time();
    numbers = _numbers;
    power = _power;
    n = 1;
    while (n * n < (int) numbers.size()) {
      n++;
    }
    assert(n * n == (int) numbers.size());
    s = *max_element(numbers.begin(), numbers.end());
    by_number.resize(s + 1);
    for (int i = 0; i < n * n; i++) {
      if (numbers[i] >= 0) {
        by_number[numbers[i]].push_back(i);
      }
    }
    start = playerRow * n + playerCol;
    nei.assign(n * n, vector<int>(4));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        for (int dir = 0; dir < 4; dir++) {
          int ti = (i + dx[dir] + n) % n;
          int tj = (j + dy[dir] + n) % n;
          nei[i * n + j][dir] = ti * n + tj;
        }
      }
    }
    vector<int> order(n * n);
    iota(order.begin(), order.end(), 0);
    vector<int> best_rm;
    vector<pair<int,int>> best_opt;
    double best_score = -2.0;
    long long ITER = 0;
    long long STEP = (int) 1e5;
    long long NEXT = STEP;
    int ITER_ID = 0;
    while (true) {
      ITER_ID++;
      ITER += STEP;
      if (ITER >= NEXT) {
        double elapsed = get_time() - start_time;
//        debug(ITER, NEXT, elapsed);
        if (elapsed > 0.9 * TL) {
          break;
        }
        NEXT += STEP;
      }
      vector<int> best_local_rm;
      double best_local_score = -2.0;
      for (int it = 0; it < 1; it++) {
        double wild_coef = rng.next_double();
        double mul_coef = 0.1 + 0.4 * rng.next_double();
        vector<int> board = numbers;
        vector<int> rm;
        int best_num = 0;
        vector<int> dfs_best_rm;
//        double utime = get_time(); //!!
        vector<int> dist(n * n, inf);
        vector<int> que;
        vector<int> inq(n * n, 0);
        int since_last = 0;
        int add_coeff = rng.next_int() % s + 1;
        function<void(int,int)> dfs = [&](int pos, int num) {
//          if (get_time() - utime > 0.02) {
//            return;
//          }
//          debug(num);
          if (num > best_num) {
            best_num = num;
            dfs_best_rm = rm;
            since_last = 0;
          } else {
            since_last++;
          }
          if (best_num == s + 1) {
            return;
          }
          if (since_last > 100) {
            return;
          }
          vector<int> pr(n * n, -1);
          fill(dist.begin(), dist.end(), inf);
          fill(inq.begin(), inq.end(), 0);
          que.assign(1, pos);
          dist[pos] = 0;
          inq[pos] = 1;
          pr[pos] = -1;
          for (int b = 0; b < (int) que.size(); b++) {
            inq[que[b]] = 0;
            if (board[que[b]] == -1 || board[que[b]] == num) {
              continue;
            }
            for (int dir = 0; dir < 4; dir++) {
              int to = nei[que[b]][dir];
              if (board[to] != 0) {
                int dt = dist[que[b]];
                if (board[to] > num) {
                  if (ITER_ID & 1) dt += s - board[to] + add_coeff; else dt += board[to];
                }
                if (dt < dist[to]) {
                  dist[to] = dt;
                  pr[to] = dir;
                  if (!inq[to]) {
                    que.push_back(to);
                    inq[to] = 1;
                  }
                }
              }
            }
          }
          vector<tuple<double,int,int>> ids;
          for (int i = 0; i < n * n; i++) {
            if (dist[i] == inf) {
              continue;
            }
            if (board[i] == num) {
              ids.emplace_back(dist[i] * (1.0 + rng.next_double() * mul_coef), rng.next_int(), i);
            }
            if (board[i] == -1) {
              ids.emplace_back(dist[i] * (1.0 + rng.next_double() * mul_coef) / wild_coef, rng.next_int(), i);
            }
          }
          if (ids.empty()) {
            return;
          }
          sort(ids.begin(), ids.end());
          for (auto &t : ids) {
            int id = get<2>(t);
            vector<pair<int,int>> board_changes;
            vector<int> cur_rm;
            int tmp = id;
            while (tmp != pos) {
              if (board[tmp] > num) {
                cur_rm.push_back(tmp);
                board_changes.emplace_back(tmp, board[tmp]);
                board[tmp] = 0;
              }
              tmp = nei[tmp][pr[tmp] ^ 2];
            }
            reverse(cur_rm.begin(), cur_rm.end());
            int old_rm_size = (int) rm.size();
            for (int x : cur_rm) {
              rm.push_back(x);
            }
            rm.push_back(id);
            if (board[id] == -1) {
              board_changes.emplace_back(id, board[id]);
              board[id] = num;
            }
            dfs(id, num + 1);
            rm.resize(old_rm_size);
            for (auto &p : board_changes) {
              board[p.first] = p.second;
            }
          }
        };
        int pos = start;
        int num = 1;
        dfs(pos, num);
        rm = dfs_best_rm;
        double score = calc(rm);
//        debug(best_num, score);
        assert(score > -0.5);
        score = best_num * 1e9 - score;
        if (score > best_local_score) {
          best_local_score = score;
          best_local_rm = rm;
        }
      }
      vector<int> rm = best_local_rm;
      double score = calc(rm);
      optimize(rm, score);
      if (score > best_score) {
        best_score = score;
        best_rm = rm;
        best_opt = opt;
      }
    }
    vector<int> rm = best_rm;
    double score = calc(rm);
    debug(score);
    string ret = build_by_rm(rm);
    debug(ret.size(), 4 * n * s); 
    if ((int) ret.size() > 4 * n * s) {
      vector<int> init_rm = rm;
      for (int i = (int) best_opt.size() - 1; i >= 0; i--) {
        init_rm.erase(init_rm.begin() + best_opt[i].first);
      }
      int low = 0, high = (int) best_opt.size();
      while (low <= high) {
        int mid = (low + high + 1) >> 1;
        rm = init_rm;
        for (int i = 0; i < mid; i++) {
          rm.insert(rm.begin() + best_opt[i].first, best_opt[i].second);
        }
        ret = build_by_rm(rm);
        if (low == high) {
          break;
        }
        if ((int) ret.size() > 4 * n * s) {
          high = mid - 1;
        } else {
          low = mid;
        }
      }
      debug(ret.size(), 4 * n * s); 
    }
    debug(get_time() - start_time);
    return ret;
  }
};
// -------8<------- end of solution submitted to the website -------8<-------

#ifdef LOCAL

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < (int) v.size(); ++i)
        cin >> v[i];
}

int main() {
    VanishingMaze vm;
    int S2;
    cin >> S2;
    vector<int> numberss(S2);
    getVector(numberss);
    int R, C;
    double P;
    cin >> R >> C >> P;

    string ret = vm.getPath(numberss, R, C, P);
    cout << ret << endl;
    cout.flush();
}

#endif
