abbrev V3 := Int × Int × Int

def signs : List (Int × Int) := [(1, 1), (1, -1), (-1, 1), (-1, -1)]

def points : List V3 :=
  ([(1, 0, 0), (0, 1, 1), (0, 1, -1), (0, 3, 2), (0, 3, -2), (1, 1, 0), (1, -1, 0)] : List V3)
  ++ (signs.map fun (e, d) => ((3 : Int), 2 * e, 2 * d))
  ++ (signs.map fun (e, d) => ((2 : Int), 3 * e, 2 * d))
  ++ (signs.map fun (e, d) => ((6 : Int), 1 * e, 2 * d))
  ++ (signs.map fun (e, d) => ((6 : Int), 5 * e, 2 * d))

def lines : List V3 :=
  ([(1, 0, 0), (0, 1, 1), (0, 1, -1), (0, 2, 3), (0, 2, -3), (1, 0, 3), (1, 0, -3)] : List V3)
  ++ (signs.map fun (e, d) => ((1 : Int), 2 * e, 2 * d))
  ++ (signs.map fun (e, d) => ((2 : Int), 6 * e, 9 * d))
  ++ (signs.map fun (e, d) => ((2 : Int), 2 * e, 1 * d))
  ++ (signs.map fun (e, d) => ((2 : Int), 2 * e, 5 * d))

def dot (p l : V3) : Int := p.1 * l.1 + p.2.1 * l.2.1 + p.2.2 * l.2.2

def onLine (p l : V3) : Bool := dot p l == 0

def cross (u v : V3) : V3 :=
  (u.2.1 * v.2.2 - u.2.2 * v.2.1, u.2.2 * v.1 - u.1 * v.2.2, u.1 * v.2.1 - u.2.1 * v.1)

def projDistinct (u v : V3) : Bool := cross u v != (0, 0, 0)

def pairs : List V3 → List (V3 × V3)
  | [] => []
  | x :: xs => xs.map (fun y => (x, y)) ++ pairs xs

def pointsPerLine (l : V3) : Nat := (points.filter (fun p => onLine p l)).length
def linesPerPoint (p : V3) : Nat := (lines.filter (fun l => onLine p l)).length
def commonLines (p q : V3) : Nat := (lines.filter (fun l => onLine p l && onLine q l)).length

def isConfig23_4 : Bool :=
  points.length == 23 && lines.length == 23
  && (pairs points).all (fun (p, q) => projDistinct p q)
  && (pairs lines).all (fun (l, m) => projDistinct l m)
  && lines.all (fun l => pointsPerLine l == 4)
  && points.all (fun p => linesPerPoint p == 4)
  && (pairs points).all (fun (p, q) => commonLines p q ≤ 1)

set_option maxRecDepth 200000 in
theorem config_23_4 : isConfig23_4 = true := by decide

#eval isConfig23_4
#print axioms config_23_4
