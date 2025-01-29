# blz78suf
blz78suf - suffix trie-based optimal LZ78 parsing for Brainfuck code
generation. Makes no assumptions on cell size, cell wrapping behaviours
or tape wrapping behaviours of the underlying interpreter. Asymptotically
optimal for highly redundant inputs.Released to the public domain by
Kamila Szewczyk - see COPYING. 

Project homepage: https://github.com/kspalaiologos/blz78suf

## Building

```
# If using a git clone (not needed for source packages), first...
$ ./bootstrap

# All...
$ ./configure
$ make
$ sudo make install
```

## Algorithm description

Standard brainfuck text generation algorithms are usually stateless (i.e.
load a byte value, output, clear, repeat) or use otherwise low order finite
state. Using [optimal constants](https://esolangs.org/wiki/Brainfuck_constants)
for such a generator renders it optimal on IID inputs. The generator is still
however unable to efficiently exploit the redundancy in such an input stemming
from the underlying probability distribution of the source.

Low-order finite state generators usually encode transitions between cell states
that are desired to be output. As an example, redundancy on the 2-wide sliding
window level can be exploited by computing a 256x256 table of optimal transition
phrases between cell values.

This however does not approximate real, variable order redundancy in the source.
Techniques that implement data compressors in brainfuck and load the compressed
data to memory, decompressing it at the runtime, generally exhibit poor performance
characteristics due to the high overhead of random memory access in Brainfuck
(fastest algorithms are quintic-time).

blz78suf operates by finding long phrases in the input that can be encoded using
procedural logic. For example:

```
void the() { printf("the"); }
int main() {
  the(); printf(" quick brown fox jumps over ");
  the(); printf(" lazy dog");
}
```

Could be shorter than:

```
int main() {
  printf("the quick brown fox jumps over the lazy dog");
}
```

Had an approperiate, sufficiently repetitive phrase been chosen. The benefit of this
approach is that we can deduplicate repeating phrases and delegate the more granular,
lower order redundancy to a stateless generator.

blz78suf builds on a stateless text generator nicked from the CodeGolf StackExchange
website for the Brainfuck Golf challenge. The same generator is used by
[copy.sh](https://copy.sh/brainfuck/text). Phrases are found by constructing the
suffix trie and ranking potential replacements by their frequency and length.
Then, the individual messages are encoded and the output is generated.

## Future improvements

blz78suf is a prototype and as such, it is not optimized for speed. The algorithm
could be sped up by using an efficient exclusion algorithm for suffix tries and
by using Ukkonen's algorithm for linear-time compressed structures. The procedural
structure of the output could be optimized by allowing phrasal chaining. Further,
a more efficient low order generator could be used. For example, such a desirable
tool would detect patterns via delta encoding and run length encoding.