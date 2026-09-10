abbrev Z17 := Int × Int
def zmul (x y : Z17) : Z17 := (x.1 * y.1 + 17 * x.2 * y.2, x.1 * y.2 + x.2 * y.1)
def zadd (x y : Z17) : Z17 := (x.1 + y.1, x.2 + y.2)
def zsub (x y : Z17) : Z17 := (x.1 - y.1, x.2 - y.2)
def zneg (x : Z17) : Z17 := (-x.1, -x.2)
def zzero : Z17 := (0, 0)
def zint (a : Int) : Z17 := (a, 0)
def zsqrt17 : Z17 := (0, 1)
def z (a b : Int) : Z17 := (a, b)

abbrev W3 := Z17 × Z17 × Z17

def signs : List (Int × Int) := [(1, 1), (1, -1), (-1, 1), (-1, -1)]
def zsc (e : Int) (x : Z17) : Z17 := (e * x.1, e * x.2)

def pointsB : List W3 :=
  ([ (z 1 0, z 0 0, z 0 0),
     (z 0 0, z 4 0, z 4 0), (z 0 0, z 4 0, z (-4) 0),                       -- (0, ±1)
     (z 0 0, z 4 0, z (-1) 1), (z 0 0, z 4 0, z 1 (-1)),                    -- (0, ±(√17-1)/4)
     (z 8 0, z 4 0, z 0 0), (z (-8) 0, z 4 0, z 0 0) ] : List W3)           -- (±2, 0)
  ++ (signs.map fun (e, d) => (zsc e (z 4 0), z 4 0, zsc d (z 3 1)))       -- (±1, ±(3+√17)/4)
  ++ (signs.map fun (e, d) => (zsc e (z 2 2), z 4 0, zsc d (z 4 0)))       -- (±(1+√17)/2, ±1)
  ++ (signs.map fun (e, d) => (zsc e (z 4 0), z 4 0, zsc d (z 5 (-1))))    -- (±1, ±(5-√17)/4)
  ++ (signs.map fun (e, d) => (zsc e (z (-2) 2), z 4 0, zsc d (z (-1) 1))) -- (±(√17-1)/2, ±(√17-1)/4)

def polar (p : W3) : W3 := (zsc 2 p.1, zsc (-4) p.2.1, zneg (zmul (z 1 1) p.2.2))

def linesB : List W3 := pointsB.map polar

def dot (p l : W3) : Z17 := zadd (zadd (zmul p.1 l.1) (zmul p.2.1 l.2.1)) (zmul p.2.2 l.2.2)
def onLine (p l : W3) : Bool := dot p l == zzero
def cross (u v : W3) : W3 :=
  (zsub (zmul u.2.1 v.2.2) (zmul u.2.2 v.2.1), zsub (zmul u.2.2 v.1) (zmul u.1 v.2.2), zsub (zmul u.1 v.2.1) (zmul u.2.1 v.1))
def isZero3 (w : W3) : Bool := w.1 == zzero && w.2.1 == zzero && w.2.2 == zzero
def projDistinct (u v : W3) : Bool := !(isZero3 (cross u v))

def pairs : List W3 → List (W3 × W3)
  | [] => []
  | x :: xs => xs.map (fun y => (x, y)) ++ pairs xs

def pointsPerLine (l : W3) : Nat := (pointsB.filter (fun p => onLine p l)).length
def linesPerPoint (p : W3) : Nat := (linesB.filter (fun l => onLine p l)).length
def commonLines (p q : W3) : Nat := (linesB.filter (fun l => onLine p l && onLine q l)).length

def isConfigB : Bool :=
  pointsB.length == 23 && linesB.length == 23
  && (pairs pointsB).all (fun (p, q) => projDistinct p q)
  && (pairs linesB).all (fun (l, m) => projDistinct l m)
  && linesB.all (fun l => pointsPerLine l == 4)
  && pointsB.all (fun p => linesPerPoint p == 4)
  && (pairs pointsB).all (fun (p, q) => commonLines p q ≤ 1)

set_option maxRecDepth 400000 in
theorem config_23_4_B : isConfigB = true := by decide

#eval isConfigB
#print axioms config_23_4_B
