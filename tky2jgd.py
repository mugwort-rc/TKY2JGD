#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

import argparse

import sys
import os

import math

import re

PAR_RE = re.compile(r"(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")

def load_parameter(par_file):
    global PAR
    PAR = {}
    for line in open(par_file):
        m = PAR_RE.match(line)
        if not m:
            continue
        PAR[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))

PAR = {}


def bilinear(lat, lon):
    """
    :type lat: float
    :type lon: float
    :rtype: float, float

    Ver.1.3  1999/3/11  (C) Mikio TOBITA 飛田幹男，国土地理院
    Bilinear interpolationをする準備とBilinear interpolationのcall
    入力座標 lat:y, lon:x   単位は度
    """
    # from Ver.1.3.78 2001/11/22
    # マイナスの緯度（南緯），経度（西経）等日本の領土外は標準地域メッシュコードが定義されていない場所では，
    # 地域毎の変換パラメータがないので，直ちにOutsideにとぶ。
    # この判定がなかったVer.1.3.77では，メッシュコードを検索に行ってしまい。見つかってしまうバグがあった。
    # なお，このバグは日本領土内では結果にまったく影響しない。
    if lat < 20 or lat > 46 or lon < 120 or lon > 154:
        return None, None

    # lat,lonからMesh Codeを作成。
    # lat,lonを含むメッシュのメッシュコードを計算
    MC0 = lat_lon2mesh_code(lat, lon)
    # MC0の東，北，東北隣のメッシュコードを計算
    MCE, MCN, MCNE = tonari_mesh_code(MC0)

    if ( MC0.mesh_code123 not in PAR or
         MCE.mesh_code123 not in PAR or
         MCN.mesh_code123 not in PAR or
         MCNE.mesh_code123 not in PAR ):
        return None, None

    dB00, dL00 = PAR[MC0.mesh_code123]
    dB10, dL10 = PAR[MCE.mesh_code123]
    dB01, dL01 = PAR[MCN.mesh_code123]
    dB11, dL11 = PAR[MCNE.mesh_code123]

    # bilinear補間。データファイルの座標系つまり日本測地系で。
    # 結果は，緯度，経度のJGD2000－日本測地系,または，Tokyo97-日本測地系
    # ModLat，ModLonは0以上1未満
    dB = interpol(dB00, dB10, dB01, dB11, MC0.mod_longitude, MC0.mod_latitude)
    dL = interpol(dL00, dL10, dL01, dL11, MC0.mod_longitude, MC0.mod_latitude)
    return dB, dL


def interpol(u1, u2, u3, u4, X0to1, Y0to1):
    """
    Ver.2.1  1999/2/4  (C) Mikio TOBITA 飛田幹男，国土地理院
    Bilinear interpolation
    X0to1, Y0to1は0以上1未満

     ^
    Y|
     u3   u4

     u1   u2  -> X

    """
    a = u1
    B = u2 - u1
    C = u3 - u1
    D = u4 - u2 - u3 + u1
    return a + B * X0to1 + C * Y0to1 + D * X0to1 * Y0to1


class MeshCode(object):
    def __init__(self, mc1, mc2, mc3, mlat=0.0, mlon=0.0):
        self.mesh_code1 = mc1
        self.mesh_code2 = mc2
        self.mesh_code3 = mc3
        self.mod_latitude = mlat
        self.mod_longitude = mlon

    @property
    def mesh_code123(self):
        return (
            self.mesh_code1 * 10000 + self.mesh_code2 * 100 + self.mesh_code3
        )

    @property
    def mesh_code_str(self):
        return "{:04d}-{:02d}-{02d}".format(
            self.mesh_code1,
            self.mesh_code2,
            self.mesh_code3
        )

def lat_lon2mesh_code(lat, lon):
    """
    :type lat: float
    :type lon: float
    :rtype: MeshCode

    Ver.2.2  2002/3/15  (C) Mikio TOBITA 飛田幹男，国土地理院
    度単位の緯度lat,経度lonを受け取り，その点が含まれるﾒｯｼｭコードを返す関数
    もし境界線上の点が与えられると，その点は点の北や東のﾒｯｼｭに属すると解釈する
    """
    # 1次メッシュコード各2桁
    lat1 = math.trunc(lat * 1.5)
    lon1 = math.trunc(lon) - 100
    # 2次メッシュコード各1桁
    lat2 = math.trunc(8.0 * (1.5 * lat - lat1))
    lon2 = math.trunc(8.0 * (lon - (lon1 + 100)))
    # 3次メッシュコード各1桁
    lat3 = math.trunc(10.0 * (12.0 * lat - 8.0 * lat1 - lat2) + 0.00000000001)
    lon3 = math.trunc(10.0 * (8.0 * (lon - (lon1 + 100)) - lon2) + 0.00000000001)
    # 上2行が微少量を加えるのは，lon=138.45のとき3次が5とならずに6となるように

    # Ver.2.1 to Ver.2.2 lat=36.0833333333333のときlat3=10になって，エラー画面が出るバグを解消。原因は微小量を加えたことの副作用。
    # 本来は微小量加算をやめるべきだが，以前の計算値との継続性を重視し，正しい緯(経)度同士の繰り上がり処理で対応する。2002/02/21
    if lat3 == 10:
        lat2 += 1
        lat3 = 0
        if lat2 == 8:
            lat1 += 1
            lat2 = 0
    if lon3 == 10:
        lon2 += 1
        lon3 = 0
        if lon2 == 8:
            lon1 += 1
            lon2 = 0

    # 3次メッシの左下（南西角）点から何度ずれているか？
    mlat = 120.0 * lat - 80.0 * lat1 - 10 * lat2 - lat3 # 余り。3次メッシの左下（南西角）点からどれくらいずれているか。0(西端)以上1(東端)未満
    mlon = 80.0 * (lon - (lon1 + 100)) - 10 * lon2 - lon3 # 余り。3次メッシの左下（南西角）点からどれくらいずれているか。0(南端)以上1(北端)未満

    # 計算値のチェック
    assert(lat2 >= 0 and lat2 <= 7)
    assert(lon2 >= 0 and lon2 <= 7)
    assert(lat3 >= 0 and lat3 <= 9)
    assert(lon3 >= 0 and lon3 <= 9)
    
    mc1 = int(lat1 * 100 + lon1)
    mc2 = int(lat2 * 10 + lon2)
    mc3 = int(lat3 * 10 + lon3)
    return MeshCode(mc1, mc2, mc3, mlat, mlon)


def tonari_mesh_code(mesh_code):
    """
    :type lat: float
    :type lon: float
    :rtype: tuple of MeshCode
    :return: (mesh_code_east, mesh_code_north, mesh_code_north_east)

    Ver.1.2  1999/2/2  (C) Mikio TOBITA 飛田幹男，国土地理院
    基本メッシュコードMeshCode1を受け取り，隣の東，北，北東のメッシュコードを返すサブルーチン
    """
    # 2 digit 1st mesh code
    lat1 = mesh_code.mesh_code1 // 100
    lon1 = mesh_code.mesh_code1 % 100
    # 1 digit 2nd mesh code
    lat2 = mesh_code.mesh_code2 // 10
    lon2 = mesh_code.mesh_code2 % 10
    # 1 digit 3rd mesh code
    lat3 = mesh_code.mesh_code3 // 10
    lon3 = mesh_code.mesh_code3 % 10

    # calc east
    lon1_tmp = int(lon1)
    lon2_tmp = int(lon2)
    lon3_tmp = int(lon3)
    if lon3_tmp != 9:
        lon3_tmp += 1
    else:
        lon3_tmp = 0
        if lon2 != 7:
            lon2_tmp += 1
        else:
            lon2_tmp = 0
            lon1_tmp += 1
    mc1 = lat1 * 100 + lon1_tmp
    mc2 = lat2 * 10 + lon2_tmp
    mc3 = lat3 * 10 + lon3_tmp
    mesh_code_east = MeshCode(mc1, mc2, mc3)

    # calc north
    lat1_tmp = int(lat1)
    lat2_tmp = int(lat2)
    lat3_tmp = int(lat3)
    if lat3_tmp != 9:
        lat3_tmp += 1
    else:
        lat3_tmp = 0
        if lat2 != 7:
            lat2_tmp += 1
        else:
            lat2_tmp = 0
            lat1_tmp += 1
    mc1 = lat1_tmp * 100 + lon1
    mc2 = lat2_tmp * 10 + lon2
    mc3 = lat3_tmp * 10 + lon3
    mesh_code_north = MeshCode(mc1, mc2, mc3)

    # north east
    mc1 = lat1_tmp * 100 + lon1_tmp
    mc2 = lat2_tmp * 10 + lon2_tmp
    mc3 = lat3_tmp * 10 + lon3_tmp
    mesh_code_north_east = MeshCode(mc1, mc2, mc3)

    return mesh_code_east, mesh_code_north, mesh_code_north_east




def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("latitude", type=float)
    parser.add_argument("longitude", type=float)
    parser.add_argument("--par", default="data/TKY2JGD.par")

    args = parser.parse_args(argv)

    if not os.path.exists(args.par):
        print("Error: Parameter file not found. \"{}\"".format(args.par), file=sys.stderr)
        return 1
    else:
        load_parameter(args.par)

    dB, dL = bilinear(args.latitude, args.longitude)
    if dB is None or dL is None:
        print(-9999.0, -9999.0)
    else:
        lon = args.longitude + dL / 3600
        lat = args.latitude + dB / 3600
        print(lat, lon)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
