// 楕円体定数（Bessel 1841）と初期化
const geo = {
  a: 6377397.155, // 長半径
  f: 1 / 299.152813, // 扁平率
  m0: 0.9999,
  Pi: Math.PI
};
geo.e2 = 2 * geo.f - geo.f ** 2; // 第一離心率²

// 系番号に対応する基準経度・緯度（北海道は12系）
function GXo(G) {
  const table = {
    12: 142.25
  };
  return table[G] || 0;
}
function GYo(G) {
  const table = {
    12: 44
  };
  return table[G] || 0;
}

// 子午線弧長 S(φ)
function S(sita) {
  const e2 = geo.e2;
  const a = geo.a;
  const terms = [
    1 + (3 / 4) * e2 + (45 / 64) * e2 ** 2,
    (3 / 4) * e2 + (15 / 16) * e2 ** 2,
    (15 / 64) * e2 ** 2,
    (35 / 512) * e2 ** 3,
    (315 / 16384) * e2 ** 4
  ];
  return a * (1 - e2) * (
    terms[0] * sita
    - terms[1] * Math.sin(2 * sita) / 2
    + terms[2] * Math.sin(4 * sita) / 4
    - terms[3] * Math.sin(6 * sita) / 6
    + terms[4] * Math.sin(8 * sita) / 8
  );
}

// φn1：垂線の足の緯度（反復計算）
function φn1(G, X) {
  let φn = GYo(G) * geo.Pi / 180;
  let Sn = S(φn);
  let M = Sn + X / geo.m0;
  let φn1;
  do {
    φn1 = φn + 2 * (Sn - M) * (1 - geo.e2 * Math.sin(φn) ** 2) ** 1.5 /
      (3 * geo.e2 * (Sn - M) * Math.sin(φn) * Math.cos(φn) * Math.sqrt(1 - geo.e2 * Math.sin(φn) ** 2)
       - 2 * geo.a * (1 - geo.e2));
    if (Math.abs(φn1 - φn) < 0.0000001) break;
    φn = φn1;
    Sn = S(φn);
  } while (true);
  return φn1;
}

// φ：緯度計算
function φ(G, X, Y) {
  const φ1 = φn1(G, X);
  const e2 = geo.e2;
  const N1 = geo.a / Math.sqrt(1 - e2 * Math.sin(φ1) ** 2);
  const t1 = Math.tan(φ1);
  const n12 = (e2 / (1 - e2)) * Math.cos(φ1) ** 2;
  const Y_ = Y / geo.m0;

  let φ = φ1
    - t1 / (2 * N1 ** 2) * (1 + n12) * Y_ ** 2
    + t1 / (24 * N1 ** 4) * (5 + 3 * t1 ** 2 + 6 * n12 - 6 * t1 ** 2 * n12 - 3 * n12 ** 2 - 9 * t1 ** 2 * n12 ** 2) * Y_ ** 4
    - t1 / (720 * N1 ** 6) * (61 + 90 * t1 ** 2 + 45 * t1 ** 4 + 107 * n12 - 162 * t1 ** 2 * n12 - 45 * t1 ** 4 * n12) * Y_ ** 6
    + t1 / (40320 * N1 ** 8) * (1385 + 3633 * t1 ** 2 + 4095 * t1 ** 4 + 1575 * t1 ** 6) * Y_ ** 8;

  return φ * 180 / geo.Pi;
}

// λ：経度計算
function λ(G, X, Y) {
  const φ1 = φn1(G, X);
  const e2 = geo.e2;
  const N1 = geo.a / Math.sqrt(1 - e2 * Math.sin(φ1) ** 2);
  const t1 = Math.tan(φ1);
  const n12 = (e2 / (1 - e2)) * Math.cos(φ1) ** 2;
  const Y_ = Y / geo.m0;

  let λ = GXo(G) * geo.Pi / 180
    + (1 / N1 / Math.cos(φ1)) * Y_
    - (1 / 6 / N1 ** 3 / Math.cos(φ1)) * (1 + 2 * t1 ** 2 + n12) * Y_ ** 3
    + (1 / 120 / N1 ** 5 / Math.cos(φ1)) * (5 + 28 * t1 ** 2 + 24 * t1 ** 4 + 6 * n12 + 8 * t1 ** 2 * n12) * Y_ ** 5
    - (1 / 5040 / N1 ** 7 / Math.cos(φ1)) * (61 + 662 * t1 ** 2 + 1320 * t1 ** 4 + 720 * t1 ** 6 * n12) * Y_ ** 7;

  return λ * 180 / geo.Pi;
}
