# TKY2JGD

国土地理院のTKY2JGDをPythonに移植したものです。

オリジナルのTKY2JGDは<http://vldb.gsi.go.jp/sokuchi/tky2jgd/download/down.cgi>からダウンロードできます。

座標変換パラメータファイル`data/TKY2JGD.par`はv2.1.1を使用しています。

## Usage

```
$ python3 tky2jgd.py 36.103774791666666 140.08785504166664
36.10696628160147 140.08457686629436
```

```
>>> import tky2jgd
>>> tky2jgd.load_parameter("data/TKY2JGD.par")
>>> lat, lon = 36.103774791666666, 140.08785504166664
>>> dB, dL = tky2jgd.bilinear(lat, lon)
>>> lat += dB / 3600
>>> lon += dL / 3600
>>> print(lat, lon)
36.10696628160147 140.08457686629436
```
