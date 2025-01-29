// ---------------------------------------------------------------------------
//      Suffix trie-based LZ78 optimal parsing and Brainfuck code
//      generation. Written on Wednesday, 29th of January 2025
//      by Kamila Szewczyk.
//
//      Makes no assumptions on cell size, cell wrapping behaviours or tape
//      wrapping behaviours of the interpreter.
//
//      See also:
//        [1] - https://copy.sh/brainfuck/text.html & improvements
//
//      TO-DO:
//        - Performance optimisations in the suffix trie search:
//          stop heuristics, better memory management, etc.
//        - Linked list-based strings for improved decomposition performance.
//        - RLE coding for uncoded phrasal sections.
//        - ...
// ---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
//      Abstract Brainfuck machine. Interprets the subset of the language
//      without I/O. Interprets the high water mark of memory usage.
//      Generates tape annihilators for arbitrary terminating programs.
// ---------------------------------------------------------------------------
class BFInt {
 public:
  std::vector<int> memory;  int ptr, mptr;
  BFInt(int cells) : memory(cells, 0), ptr(0), mptr(0) {}
  void simulate(const std::string & program) {
    for (int i = 0; i < program.size(); i++) {
      char c = program[i];
      switch (c) {
        case '>': ptr++; break; case '<': ptr--; break;
        case '+': memory[ptr]++; break;
        case '-': memory[ptr]--; break;
        case '[': {
          if (!memory[ptr]) {
            for (int depth = 1; depth != 0; ) {
              i++;
              if (program[i] == '[') depth++;
              else if (program[i] == ']') depth--;
            }
          }
          break;
        }
        case ']': {
          if (memory[ptr]) {
            for (int depth = 1; depth != 0; ) {
              i--;
              if (program[i] == '[') depth--;
              else if (program[i] == ']') depth++;
            }
          }
          break;
        }
      }
      if (ptr > mptr) mptr = ptr;
    }
  }
  std::string craftAnnihilator() {
    std::string annihilator = "";
    annihilator.reserve(32);
    while (true) {
      switch (memory[ptr]) {
        case 0: break;
        case 1: annihilator += "-"; break;
        case 2: annihilator += "--"; break;
        default: annihilator += "[-]"; break;
      }
      memory[ptr] = 0;
      int i;
      for (i = ptr; memory[i] == 0 && i < memory.size(); i++);
      if (i != memory.size()) {
        annihilator += std::string(i - ptr, '>');
        ptr = i; continue;
      }
      for (i = ptr; memory[i] == 0 && i >= 0; i--);
      if (i != -1) {
        annihilator += std::string(ptr - i, '<');
        ptr = i; continue;
      }
      annihilator += std::string(ptr, '<');
      ptr = 0; return annihilator;
    }
  }
};

// ---------------------------------------------------------------------------
//      Non-quadratic approximator for uncoded phrase generation.
// ---------------------------------------------------------------------------
class BFGenApprox {
 private:
  int G[256][256];
 public:
  BFGenApprox() {
    for (int x = 0; x < 256; x++) {
      for (int y = 0; y < 256; y++) {
        int delta = y - x;
        if (delta > 128) delta -= 256;
        if (delta < -128) delta += 256;
        G[x][y] = delta >= 0 ? delta : -delta;
      }
    }
    bool iter = true;
    while (iter) {
      iter = false;
      for (int x = 0; x < 256; x++) {
        for (int n = 1; n < 40; n++) {
          for (int d = 1; d < 40; d++) {
            int j = x; int y = 0;
            for (int i = 0; i < 256; i++) {
              if (j == 0) break;
              j = (j - d + 256) % 256;
              y = (y + n) % 256;
            }
            if (j == 0) {
              int s = 5 + d + n;
              if (s < G[x][y]) {
                G[x][y] = s;
                iter = true;
              }
            }
            j = x; y = 0;
            for (int i = 0; i < 256; i++) {
              if (j == 0) break;
              j = (j + d) % 256;
              y = (y - n + 256) % 256;
            }
            if (j == 0) {
              int s = 5 + d + n;
              if (s < G[x][y]) {
                G[x][y] = s;
                iter = true;
              }
            }
          }
        }
      }
      for (int x = 0; x < 256; x++) {
        for (int y = 0; y < 256; y++) {
          for (int z = 0; z < 256; z++) {
            if (G[x][z] + G[z][y] < G[x][y]) {
              G[x][y] = G[x][z] + G[z][y];
              iter = true;
            }
          }
        }
      }
    }
  }
  size_t phrase_len(const std::string & s) {
    int lastc = 0, gen = 0;
    for (char c : s) {
      int a = G[lastc][c], b = G[0][c];
      if (a + 3 <= b) gen += a + 1;
      else            gen += b + 4;
      lastc = c;
    }
    gen += 3; return gen;
  }
};

// ---------------------------------------------------------------------------
//      Low-level, precise text generation for uncoded phrases.
// ---------------------------------------------------------------------------
enum class Direction { LEFT, RIGHT };
Direction opposite(Direction d) {
  return d == Direction::LEFT ? Direction::RIGHT : Direction::LEFT;
}
class Transition {
 public:
  std::string code; Direction startD, endD;
  Transition() : code(""), startD(Direction::LEFT), endD(Direction::LEFT) {}
  Transition(std::string code, Direction startD, Direction endD)
    : code(code), startD(startD), endD(endD) {}
  Transition(const Transition & t)
    : code(t.code), startD(t.startD), endD(t.endD) {}
  Transition(Transition && t)
    : code(std::move(t.code)), startD(t.startD), endD(t.endD) {}
  Transition & operator=(const Transition & t) {
    code = t.code; startD = t.startD; endD = t.endD; return *this;
  }
  Transition plus(Transition t) const {
    if (endD == t.startD) {
      return Transition(std::move(code + t.code), t.startD, t.endD);
    } else {
      Direction tOriginalStartD = t.startD; t.reverse_in_place();
      return Transition(std::move(code + t.code), tOriginalStartD, t.endD);
    }
  }
  void reverse_in_place() {
    for (int i = 0; i < code.size(); i++) {
      if (code[i] == '>')      code[i] = '<';
      else if (code[i] == '<') code[i] = '>';
    }
    startD = opposite(startD); endD = opposite(endD);
  }
  int size() const { return code.size(); }
};
class BFGen {
 private:
  Transition list[256][256];
  std::string generateFromCache(const std::string & s,
      const std::vector<bool> & caches) {
    std::vector<Transition> trans;
    trans.reserve(s.size());
    Direction d = Direction::LEFT;
    char last = 0;
    std::vector<char> cache;
    for (int index = 0; index < s.size(); index++) {
      char c = s[index];  Transition t = list[last][c];
      if (t.startD != d) t.reverse_in_place();
      int ci = std::find(cache.begin(), cache.end(), c) - cache.begin();
      int ims = (cache.size() - ci) * 2 + (d == Direction::RIGHT ? 2 : 0);
      if (ci != cache.size() && ims <= t.size())
        t = Transition(std::to_string(cache.size() - ci), t.startD, t.startD);
      else {
        t = Transition(t.code + ".", t.startD, t.endD); d = t.endD; last = c;
      }
      if (caches[index]) {
        cache.push_back(c); d = Direction::LEFT; last = 0;
      }
      trans.emplace_back(std::move(t));
    }
    bool reverse = false;
    for (int i = trans.size() - 1; i >= 0; i--) {
      if (caches[i]) reverse = trans[i].endD == Direction::RIGHT;
      if (reverse) trans[i].reverse_in_place();
    }
    if (trans.empty()) return "";
    std::string code = "";
    for (int i = 0; i < trans.size(); i++) {
      Transition t = trans[i];
      if ((i == 0 && t.startD != Direction::LEFT)
      || (i != 0 && trans[i - 1].endD != t.startD))
        code += ">";
      if (std::isdigit(t.code[0])) {
        int num = std::stoi(t.code) + (t.startD == Direction::RIGHT);
        code += std::string(num, '<') + "." + std::string(num, '>');
      } else
        code += t.code;
      if (caches[i]) code += ">";
    }
    while (true) {
      std::string newCode = code; size_t pos = 0;
      while ((pos = newCode.find("><", pos)) != std::string::npos)
        newCode.replace(pos, 2, "");
      if (code == newCode) break; code = newCode;
    }
    return code;
  }
  std::string generate_internal(const std::string & s, int cells) {
    std::vector<bool> caches(s.size(), false);
    std::string currentCost = generateFromCache(s, caches);
    std::map<char, int> ch;
    if (cells > 2) {
      for (int i = 0; i < s.size(); i++) {
        char c = s[i]; std::vector<bool> nc = caches;
        for (int j = 0; j < s.size(); j++) {
          if (s[j] == c) nc[j] = false;
        }
        char lowest = 0;
        if (std::count(nc.begin(), nc.end(), true) >= cells - 2) {
          lowest = std::min_element(ch.begin(), ch.end(),
              [](auto & a, auto & b) { return a.second < b.second; }
            )->first;
          auto iter = std::find(nc.begin(), nc.end(), true);
          nc[std::distance(nc.begin(), iter)] = false;
        }
        nc[i] = true;
        std::string newCost = generateFromCache(s, nc);
        if (newCost.size() < currentCost.size()) {
          ch[c] = currentCost.size() - newCost.size();
          currentCost = newCost; caches = nc;
          if (lowest != 0) ch.erase(lowest);
        }
      }
    }
    return currentCost;
  }
  int grade(int n, int base) {
    int sp = 0, norm = 0;
    while (n > 0) {
      sp++; norm += n % base; n = n / base;
    }
    return norm + (6 + base) * sp + ((sp % 2 == 1) ? 4 : 0);
  }
  int best_base(int n) {
    int v = 0, b = 0;
    for (int i = 2; i <= 60; i++) {
      int cv = grade(n, i);
      if (v == 0 || v > cv) {
        v = cv; b = i;
      }
    }
    return b;
  }
 public:
  BFGen() {
    for (int x = 0; x < 256; x++) {
      for (int y = 0; y < 256; y++) {
        int delta = y - x;
        std::string code = "";
        if (delta > 0)
          code = std::string(delta, '+');
        else if (delta < 0)
          code = std::string(-delta, '-');
        list[x][y] = Transition(code, Direction::LEFT, Direction::LEFT);
      }
    }
    for (int x = 0; x < 256; x++) {
      for (int n = 1; n <= 39; n++) {
        for (int d = 1; d <= 39; d++) {
          int j = x; int y = 0;
          for (int i = 0; i < 256; i++) {
            if (j == 0 || j - d < 0 || y + n > 255) break;
            j = (j - d + 256) & 255;
            y = (y + n) & 255;
          }
          if (j == 0) {
            std::string s =
              "[" + std::string(d, '-') + ">" + std::string(n, '+') + "<]>";
            if (s.size() < list[x][y].size())
              list[x][y] = Transition(s, Direction::LEFT, Direction::RIGHT);
          }
          j = x; y = 0;
          for (int i = 0; i < 256; i++) {
            if (j == 0 || y - n < 0 || j + d > 255) break;
            j = (j + d) & 255;
            y = (y - n + 256) & 255;
          }
          if (j == 0) {
            std::string s =
              "[" + std::string(d, '+') + ">" + std::string(n, '-') + "<]>";
            if (s.size() < list[x][y].size())
              list[x][y] = Transition(s, Direction::LEFT, Direction::RIGHT);
          }
        }
      }
    }
    for (int i = 0; i < 256; i++)
      if (list[i][0].size() > 3)
        list[i][0] = Transition("[-]", Direction::LEFT, Direction::LEFT);
    bool change = true;
    while (change) {
      change = false;
      for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
          for (int k = 0; k < 256; k++) {
            if (list[i][j].size() + list[j][k].size() < list[i][k].size()) {
              list[i][k] = list[i][j].plus(list[j][k]); change = true;
            }
          }
        }
      }
    }
  }
  std::string generate(const std::string & s, int cells_max = 16) {
    std::string code = generate_internal(s, cells_max);
    BFInt bf(cells_max); bf.simulate(code);
    cells_max = bf.mptr + 1;
    std::string shortest_code = code + bf.craftAnnihilator();
    while (cells_max > 2) {
      std::string new_code = generate_internal(s, cells_max - 1);
      BFInt bf(cells_max - 1); bf.simulate(new_code);
      new_code += bf.craftAnnihilator();
      if (new_code.size() < shortest_code.size())
        shortest_code = new_code;
      cells_max--;
    }
    return shortest_code;
  }
  std::string gen_constant(int n) {
    std::vector<int> stack;
    std::string out = ">";
    int flip = 1;
    if (n < 12)
      return std::string(n, '+');
    int base = best_base(n);
    while (n > 0) {
      stack.push_back(n % base);
      n = n / base;
    }
    while (!stack.empty()) {
      int top = stack.back();
      stack.pop_back();
      int bc = base;
      out += std::string(top, '+');
      if (!stack.empty()) {
        if (!flip)
          out += "[>" + std::string(bc, '+') + "<-]>";
        else
          out += "[<" + std::string(bc, '+') + ">-]<";
      }
      flip = !flip;
    }
    if (!flip) out += "[-<+>]<";
    return out;
  }
};

// ---------------------------------------------------------------------------
//      Suffix trie for ranked LZ78 parsing.
// ---------------------------------------------------------------------------
class SuffixTrie {
 private:
  class SuffixTrieNode {
   public:
    std::unordered_map<char, SuffixTrieNode *> children;
    int counter, cache;
    SuffixTrieNode() : counter(0), cache(0) {}
  };
  SuffixTrieNode * root;
  void buildTrie(const std::string & text) {
    int n = text.length();
    for (int i = 0; i < n; i++) {
      SuffixTrieNode * cn = root;
      for (int j = i; j < n; j++) {
        char ch = text[j];
        if (cn->children.find(ch) == cn->children.end()) {
          cn->children[ch] = new SuffixTrieNode();
        }
        cn = cn->children[ch]; cn->counter++;
      }
    }
  }
  void annotateCounters(SuffixTrieNode * node) {
    if (node->children.empty()) {
      node->cache = 1; return;
    }
    for (auto & pair : node->children) {
      annotateCounters(pair.second);
      node->cache += pair.second->cache;
    }
  }
 public:
  SuffixTrie(const std::string & text) {
    this->root = new SuffixTrieNode();
    buildTrie(text); annotateCounters(root);
  }
  ~SuffixTrie() { deleteTrie(root); }
  template <typename F> std::pair<std::string, double> findMaxString(F rate,
      double constant, double constant2) {
    double maxValue = -std::numeric_limits<double>::infinity();
    std::string bestString;
    auto dfs = [&](this auto const & dfs, SuffixTrieNode * node,
        std::string currentString) -> void {
      if (!currentString.empty() && node->counter > 1
      && currentString.find('\x01') == std::string::npos) {
        int countX = node->cache; double rateX = rate(currentString);
        double value = countX * rateX - countX * constant2 - constant;
        if (value > maxValue) {
          maxValue = value; bestString = currentString;
        }
      }
      for (auto & pair : node->children)
        dfs(pair.second, currentString + pair.first);
    };
    dfs(root, ""); return { bestString, maxValue };
  }
 private:
  void deleteTrie(SuffixTrieNode * node) {
    for (auto & pair : node->children)
      deleteTrie(pair.second);
    delete node;
  }
};

// ---------------------------------------------------------------------------
//      Code generation and procedural system.
// ---------------------------------------------------------------------------
class CodeGen {
 private:
  std::string replacement_program(BFGen & gen, const std::string & msg) {
    return "-[-<+>>>>+<<<]<[->+<]>>>>>[-]<[>+<[-]]+>[<->-]<[-" +
           gen.generate(msg) + "<<[->+<]>>]<<<";
  }
 public:
  std::string gen(BFGen & gen, const std::vector<std::string> & replacements,
      const std::vector<std::string> & chunks) {
    int ret = 0; std::string new_code;
    new_code = ">" + gen.gen_constant(replacements.size() + 1) + "[";
    for (const auto & el : replacements)
      new_code += replacement_program(gen, el);
    new_code += "-[-<+>>>>+<<<]<[->+<]>>>>>[-]<[>+<[-]]+>[<->-]<[-";
    ret++;
    for (const auto & el : chunks) {
      auto iter = std::find(replacements.begin(), replacements.end(), el);
      if (iter != replacements.end()) {
        int c1 = (ret++) + replacements.size() + 1;
        int c2 = std::distance(replacements.begin(), iter) + 1;
        new_code += "<<" + gen.gen_constant(c1) + ">" + gen.gen_constant(c2) +
                    ">]<<<-[-<+>>>>+<<<]<[->+<]>>>>>[-]<[>+<[-]]+>[<->-]<[-";
      } else new_code += gen.generate(el);
    }
    return new_code + "]<<<[-]>>[-<<+>>]<<]";
  }
};

// ---------------------------------------------------------------------------
//      Constrained, optimal phrasal LZ78 parsing of text.
// ---------------------------------------------------------------------------
template <typename F>
std::vector<std::string> parse(std::string text, F score, int max,
    int c1, int c2) {
  if (text.contains("\x01"))
    throw std::runtime_error("Text contains forbidden characters.");
  std::vector<std::string> replacements;
  double maxValue = std::numeric_limits<double>::infinity();
  while (maxValue > 0 && replacements.size() < max) {
    SuffixTrie suffixTrie(text); std::string bestString;
    std::tie(bestString, maxValue) = suffixTrie.findMaxString(score, c1, c2);
    if (maxValue > 0) replacements.push_back(bestString);
    for (size_t pos = 0;
        (pos = text.find(bestString, pos)) != std::string::npos;
        pos++) text.replace(pos, bestString.length(), "\x01");
  }
  return replacements;
}
void splitText(std::vector<std::string> & acc, const std::string & text,
    const std::vector<std::string> & words) { // TO-DO: `words' as a map.
  if (text.empty()) return;
  std::string longest = "";
  for (const auto & word : words)
    if (text.rfind(word, 0) == 0 && word.length() > longest.length())
      longest = word;
  if (longest.empty()) {
    if (acc.empty() || find(words.begin(), words.end(), acc.back()) != words.end())
      acc.push_back(std::string(1, text[0]));
    else acc.back() += text[0];
    return splitText(acc, text.substr(1), words);
  } else {
    acc.push_back(longest);
    return splitText(acc, text.substr(longest.length()), words);
  }
}

// ---------------------------------------------------------------------------
//      Command-line stub.
// ---------------------------------------------------------------------------
int main(int argc, char * argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: blz78suf <text_file> <c1> <c2> <max>" << std::endl;
    return 1;
  }
  std::ifstream file(argv[1]);
  if (!file.is_open()) {
    std::cerr << "Error: could not open file " << argv[1] << std::endl;
    return 1;
  }
  int c1 = std::stoi(argv[2]), c2 = std::stoi(argv[3]),
      max = std::stoi(argv[4]);
  std::string text = std::string(std::istreambuf_iterator<char>(file),
                                 std::istreambuf_iterator<char>());
  BFGen bfgen; std::string naive_gen = bfgen.generate(text);
  std::cout << "Naive (" << naive_gen.length() << " bytes): "
            << naive_gen << std::endl;
  BFGenApprox approx;
  std::vector<std::string> replacements =
      parse(text, [&](const std::string & s) { return approx.phrase_len(s); },
            c1, c2, max);
  std::vector<std::string> acc;
  splitText(acc, text, replacements);
  std::cout << "LZ78 parsing produced " << replacements.size()
            << " phrases and " << acc.size() << " tokens." << std::endl;
  CodeGen cg; std::string clever_gen = cg.gen(bfgen, replacements, acc);
  std::cout << "Clever (" << clever_gen.length() << " bytes): "
            << clever_gen << std::endl;
}
