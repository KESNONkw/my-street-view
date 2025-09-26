// Entry point: 電柱番号 → 緯度・経度
function getLatLon(denchuNumber) {
  const { X, Y } = getXY(denchuNumber);
  const zone = 12; // 北海道12系
  return {
    lat: calcLatitude(zone, X, Y),
    lng: calcLongitude(zone, X, Y)
  };
}

// 楕円体定数（Bessel 1841）
const GEO = {
  a: 6377397.155,
  f: 1 / 299.152813,
  m0: 0.9999,
  Pi: Math.PI
};
GEO.e2 = 2 * GEO.f - GEO.f ** 2;

// 系番号 → 基準経度・緯度
function getCentralLongitude(G) {
  return { 12: 142.25 }[G] || 0;
}
function getCentralLatitude(G) {
  return { 12: 44 }[G] || 0;
}

// 子午線弧長 S(φ)
function calcMeridianArc(phi) {
  const { a, e2 } = GEO;
  const terms = [
    1 + (3 / 4) * e2 + (45 / 64) * e2 ** 2,
    (3 / 4) * e2 + (15 / 16) * e2 ** 2,
    (15 / 64) * e2 ** 2,
    (35 / 512) * e2 ** 3,
    (315 / 16384) * e2 ** 4
  ];
  return a * (1 - e2) * (
    terms[0] * phi
    - terms[1] * Math.sin(2 * phi) / 2
    + terms[2] * Math.sin(4 * phi) / 4
    - terms[3] * Math.sin(6 * phi) / 6
    + terms[4] * Math.sin(8 * phi) / 8
  );
}

// 垂線の足 φn1（反復）
function calcFootLatitude(G, X) {
  let phi = getCentralLatitude(G) * GEO.Pi / 180;
  let S = calcMeridianArc(phi);
  const M = S + X / GEO.m0;
  let nextPhi;
  do {
    nextPhi = phi + 2 * (S - M) * (1 - GEO.e2 * Math.sin(phi) ** 2) ** 1.5 /
      (3 * GEO.e2 * (S - M) * Math.sin(phi) * Math.cos(phi) * Math.sqrt(1 - GEO.e2 * Math.sin(phi) ** 2)
       - 2 * GEO.a * (1 - GEO.e2));
    if (Math.abs(nextPhi - phi) < 1e-7) break;
    phi = nextPhi;
    S = calcMeridianArc(phi);
  } while (true);
  return nextPhi;
}

// 緯度 φ
function calcLatitude(G, X, Y) {
  const phi1 = calcFootLatitude(G, X);
  const { a, e2, m0, Pi } = GEO;
  const N1 = a / Math.sqrt(1 - e2 * Math.sin(phi1) ** 2);
  const t1 = Math.tan(phi1);
  const n12 = (e2 / (1 - e2)) * Math.cos(phi1) ** 2;
  const Y_ = Y / m0;

  const phi = phi1
    - t1 / (2 * N1 ** 2) * (1 + n12) * Y_ ** 2
    + t1 / (24 * N1 ** 4) * (5 + 3 * t1 ** 2 + 6 * n12 - 6 * t1 ** 2 * n12 - 3 * n12 ** 2 - 9 * t1 ** 2 * n12 ** 2) * Y_ ** 4
    - t1 / (720 * N1 ** 6) * (61 + 90 * t1 ** 2 + 45 * t1 ** 4 + 107 * n12 - 162 * t1 ** 2 * n12 - 45 * t1 ** 4 * n12) * Y_ ** 6
    + t1 / (40320 * N1 ** 8) * (1385 + 3633 * t1 ** 2 + 4095 * t1 ** 4 + 1575 * t1 ** 6) * Y_ ** 8;

  return phi * 180 / Pi;
}

// 経度 λ
function calcLongitude(G, X, Y) {
  const phi1 = calcFootLatitude(G, X);
  const { a, e2, m0, Pi } = GEO;
  const N1 = a / Math.sqrt(1 - e2 * Math.sin(phi1) ** 2);
  const t1 = Math.tan(phi1);
  const n12 = (e2 / (1 - e2)) * Math.cos(phi1) ** 2;
  const Y_ = Y / m0;

  const lambda = getCentralLongitude(G) * Pi / 180
    + (1 / N1 / Math.cos(phi1)) * Y_
    - (1 / 6 / N1 ** 3 / Math.cos(phi1)) * (1 + 2 * t1 ** 2 + n12) * Y_ ** 3
    + (1 / 120 / N1 ** 5 / Math.cos(phi1)) * (5 + 28 * t1 ** 2 + 24 * t1 ** 4 + 6 * n12 + 8 * t1 ** 2 * n12) * Y_ ** 5
    - (1 / 5040 / N1 ** 7 / Math.cos(phi1)) * (61 + 662 * t1 ** 2 + 1320 * t1 ** 4 + 720 * t1 ** 6 * n12) * Y_ ** 7;

  return lambda * 180 / Pi;
}

// 電柱番号を分解
function parseDenchuNumber(number) {
  return {
    kuY: parseInt(number.slice(2, 3)),
    kuX: parseInt(number.slice(3, 4)),
    zuY: parseInt(number.slice(4, 5)),
    zuX: parseInt(number.slice(5, 6)),
    banY: parseInt(number.slice(6, 7)),
    banX: parseInt(number.slice(7, 8)),
    no: parseInt(number.slice(8, 10)),
    gou: parseInt(number.slice(10, 12))
  };
}

// 各階層のサイズ（画42基準）
const GRID_SIZE = {
  kuY: 9373.875,
  kuX: 10999.125,
  zuY: 937.3875,
  zuX: 1099.9125,
  banY: 93.73875,
  banX: 109.99125
};

// 電柱番号 → X,Y距離
function getXY(denchuNumber) {
  const d = parseDenchuNumber(denchuNumber);
  const X = d.kuY * GRID_SIZE.kuY + d.zuY * GRID_SIZE.zuY + d.banY * GRID_SIZE.banY + d.no;
  const Y = d.kuX * GRID_SIZE.kuX + d.zuX * GRID_SIZE.zuX + d.banX * GRID_SIZE.banX + d.gou;
  return { X, Y };
}
