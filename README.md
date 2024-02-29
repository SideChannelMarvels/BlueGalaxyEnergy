# Blue Galaxy Energy

*Hologram: Shut up! You do not know the power of the Blue Galaxy Energy! Also known as the "B.G.E"*
*Mr. Whereabout: The Loss, Part III, Volume I*

BlueGalaxyEnergy is a tool to perform the so-called BGE attack described in

- *Cryptanalysis of a White Box AES Implementation*, Olivier Billet, Henri Gilbert, Charaf Ech-Chatbi

with the optimizations proposed in:

- *Improved cryptanalysis of an AES implementation*, Ludo Tolhuizen
- *Revisiting the BGE Attack on a White-Box AES Implementation*, Yoni De Mulder, Peter Roelse, Bart Preneel

## Compile

To compile and install the project, install gmp and ntl libraries and development headers (available in your OS package manager), e.g.

```bash
$ sudo apt install libgmp-dev libntl-dev
```
or
```bash
$ sudo pacman -S gmp ntl
```

then compile and install the package locally:

```bash
$ pip install bluegalaxyenergy
```

## Test

```bash
$ python3 -m bluegalaxyenergy --selftest
```

## Run the attack

```python
from bluegalaxyenergy import WhiteBoxedAES, BGE

class MyWhiteBoxedAES(WhiteBoxedAES):

    def __init__(self, ...):
        # TODO

    def isEncrypt(self):
        # return True if the whitebox is an encryption whitebox, False otherwise
        return True

    def getRoundNumber(self):
        # return the number of rounds of the whitebox (10 for AES128,
        #   12 for AES192 and 14 for AES256)
        return 10

    def applyRound(self, data, roundN):
        # Apply a round of the whitebox on a buffer
        # [param] data    a buffer of 16 bytes (type bytes)
        # [param] roundN  the round number to apply (int in the range [0, self.getRoundNumber()) )
        # return  16 bytes of the encrypted data by the round
        return ... # TODO

mywb = MyWhiteBoxedAES(...)

# run the attack
bge = BGE(mywb)
bge.run()

# extract the key from the available roundKey
key = bge.computeKey()
if key is not None:
    print("key:", key.hex())
```

By default, the method `run` extracts all rounds except for the last one, which
lacks a MixColumns operation. You can restrict the rounds used for the attack
using the `roundList` option: `run(roundList = [4, 5, 6, 7, 8])`.
However, the attack requires a minimum of three consecutive rounds to extract a single round key.

- AES128 necessitates a minimum of three consecutive rounds for key extraction.
- AES192 and AES256 require a minimum of four consecutive rounds.

If the whitebox's intermediary state is shuffled, specify `shuffle=True` in the `run` method.
When shuffled, the number of required rounds increases to four for AES128 and AES192, and five for AES256.
Insufficient rounds will result in no key being found or a dictionary containing 16 possible keys.

The bge class can additionally extract the encoding of the intermediary round (with `bge.getEncoding()`)
and the permutation of the intermediary round if the whitebox is shuffled (with `bge.getShuffle()`).
For more information on these elements, refer to the test implemented in [`__main__.py`](src/bluegalaxyenergy/__main__.py).

## About

### Authors and Contributors

Initial Authors and Contributors:

- Laurent Grémy
- Nicolas Surbayrole
- Philippe Teuwen

For subsequent contributions, see the git projet history.

### Copyright

[Quarkslab](https://www.quarkslab.com)

### License

BlueGalaxyEnergy is provided under the [Apache 2.0 license](LICENSE.txt).
