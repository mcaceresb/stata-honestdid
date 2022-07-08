HonestDiD
=========

HonestDiD translation to Stata xx

`version 0.1.0 07Jul2022` | [Installation](#installation) | [Usage](#usage) | [Examples](#examples)

### Installation

From Stata

```stata
local github "https://raw.githubusercontent.com"
cap noi net uninstall honestdid
net install honestdid, from(`github'/mcaceresb/stata-honestdid/main/)
```

You can also clone or download the code manually, e.g. to
`stata-honestdid-main`, and install from a local folder:

```stata
cap noi net uninstall honestdid
net install honestdid, from(`c(pwd)'/stata-honestdid-main)
```

### Usage

```
honestdid
honestdid, [b() vcov() coefplot ...]
help honestdid
```

### Examples

See `./test/test-replication.do`; compare to `./test/test-replication.R`.
