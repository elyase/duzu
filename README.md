# duzu

Python port of the 10 parameter [Duflo and Zucker](http://amdc.in2p3.fr/web/dz.html) (J. Duflo and A.P. Zuker
*Phys. Rev. C* **52**, 1995) formula.


## Usage:


```python
>>> from duzu import duzu10
>>> duzu10(82,126)         # Calculates the mass excess for (N,Z) = (82,126)
-22.583689938619045
```

## Install

Download the script and place it somewhere your program can see it for example:

```bash
curl -O https://raw.github.com/elyase/duzu/master/duzu.py
```

## Requirements

* numpy
* pytest (optional for testing)

## Testing

Make sure you have py.test installed then run:

```bash
py.test
```

## License

The MIT License (MIT)

Copyright (c) 2013, Yaser Martinez Palenzuela

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

