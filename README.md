# DitheringScript

## 重要なお知らせ

**プログラムを一新したバージョンを公開したので、こちらのリポジトリはアーカイブにしました**

**新しいバージョンは[コチラから](https://github.com/Totti95U/aviutl-dither-script)**

---

配列ディザリング (パターンディザリング) を行う AviUtl のスクリプトです。

ダウンロードは[コチラ](https://github.com/Totti95U/DitheringScript/releases)

現在用意している効果は以下の通りです。

---

## @Simple Dithering
白、黒、赤、青、緑、シアン、マゼンダ、イエローの8色によるディザ処理を行います。

![simple dithering](https://github.com/Totti95U/DitheringScript/blob/images/Simple_Dithering.png)

---
## @DitherComposer
透明度に対して配列ディザ処理を行います。(Adobe Illustrater の ディザ合成と同等)

![dither composer](https://github.com/Totti95U/DitheringScript/blob/images/Dither_Composer.png)

---
## @Color Degreasing
簡易な中央値分割法によりオブジェクトの色数を減らします。

![color degreasing](https://github.com/Totti95U/DitheringScript/blob/images/Color_Degreasing.png)

---
## @Color Degreasing to W3C pallet
オブジェクトをW3C標準 基本16色に減色します。

![color degreasing to @3C pallet](https://github.com/Totti95U/DitheringScript/blob/images/Color_Degreasing_to_W3C.png)

---
# 調整項目について
各効果にある調整項目について説明します。

## ドットサイズ
画像の画素数を調整します。(0：最小　100：最大)

## 色数
いくつの色に減色するかを選択します。
