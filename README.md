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

By default, the method `run()` will extract all the rounds (except the last one who
doesn't have a mixColumns). You can limit which round to use for the attack
with the option `roundList` : `bge.run(roundList = [4,5,6,7,8])`. However, the
attack needs three consecutive rounds in order to extract one round key.
For AES128, a minimum of three consecutive rounds is needed to extract the key.
For AES192 and AES256, the minimum is four consecutive rounds.

## Limitations

The implementation can be used on AES-whiteboxes which encrypt, not on those
which decrypt.

The implementation does not cover the randomization in the order of the bytes
of the intermediate results in AES, mentioned in De Mulder paper.

## About

### Authors and Contributors

Initial Authors and Contributors:

- Laurent Grémy
- Nicolas Surbayrole
- Philippe Teuwen

For next contributions, see the git projet history.

### Copyright

[Quarkslab](https://www.quarkslab.com)

### License

BlueGalaxyEnergy is provided under the [Apache 2.0 license](LICENSE.txt).
