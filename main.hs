import Codec.Picture

data Scene = Mirror    Scene Plane
           | Reflect   Scene Plane
           | Translate Scene Vector
           | Rotate    Scene Quaternion
           | ScaleO    Scene Double
           | Triangles [Triangle]

type Vector = (Double, Double, Double)
type Quaternion = (Double, Vector)
type Triangle = (Vector, Vector, Vector)
type Plane = (Vector, Vector)
type Vect2d = (Double, Double)
type Matrix = (Vect2d, Vect2d)
type CompressedScene = ([Triangle], [Plane])     --includes a set of triangles and a set of mirror planes
type Line = (Vect2d, Vect2d)

--Warning: there are a lot of unused functions and very few comments.

matMult :: Matrix -> Matrix -> Matrix
matMult ((a11,a21),(a12,a22)) ((b11,b21),(b12,b22)) = ((a11*b11 + a12*b21, a21*b11 + a22*b21), (a11*b12 + a12*b22, a21*b12 + a22*b22))

matVecMult :: Matrix -> Vect2d -> Vect2d
matVecMult (c1,c2) (v1,v2) = (v1 *.. c1) +.. (v2 *.. c2)

matInverse :: Matrix -> Matrix
matInverse ((a,b),(c,d)) = ( (d / det, -b / det), (-c / det, a / det) )
  where det = (a*d) - (b*c)

cross :: Vector -> Vector -> Vector
cross (x1, y1, z1) (x2, y2, z2) = (y1*z2 - y2*z1, z1*x2 - z2*x1, x1*y2 - x2*y1)

dot :: Vector -> Vector -> Double
dot (x1, y1, z1) (x2, y2, z2) = x1*x2 + y1*y2 + z1*z2

(*.) :: Double -> Vector -> Vector
s *. (x, y, z) = (s*x, s*y, s*z)

(+.) :: Vector -> Vector -> Vector
(x1,y1,z1) +. (x2,y2,z2) = (x1+x2, y1+y2, z1+z2)

project :: Vector -> Vector -> Vector
project v1 v2 = ( (dot v1 v2) / (dot v1 v1) ) *. v1

vLength :: Vector -> Double
vLength v = sqrt (dot v v)

qMult :: Quaternion -> Quaternion -> Quaternion
qMult (r1, v1) (r2, v2) = (r1*r2 - (dot v1 v2), (r1 *. v2) +. (r2 *. v1) +. (cross v1 v2))

constructRotation :: Vector -> Double -> Quaternion
constructRotation unitAxis angle = (cos (angle/2), (sin (angle/2)) *. unitAxis)

inv :: Quaternion -> Quaternion
inv (s, v) = qMult (s, (-1) *. v) (1 / (s*s + dot v v), (0.0, 0.0, 0.0))

dmod :: Double -> Double -> Double
dmod x y = x - (y * (fromIntegral ( truncate (x/y + 0.0000000001))))    -- The small addition is to ensure that with x/y = integer, it never goes an integer lower than it should.

rotateVect :: Quaternion -> Vector -> Vector
rotateVect q v = (\(scal,vect) -> vect) ( qMult (qMult q (0,v)) (inv q) )

factorToMeetPlane :: Plane -> Vector -> Double
factorToMeetPlane (plPoint, vPerp) vect = (dot plPoint vPerp) / (dot vect vPerp)

planeReflect :: Plane -> Vector -> Vector
planeReflect (plPoint, vPerp) point = point +. ( (-2) *. ( project vPerp (point +. ((-1) *. plPoint)) ) )      -- point - 2*(  project vPerp (plPoint to point)  )

planePlaneReflect :: Plane -> Plane -> Plane
planePlaneReflect plane (plPoint, vPerp) = (planeReflect plane plPoint, planeReflect ((\(plPoint,vPerp) -> ((0,0,0), vPerp)) plane) vPerp)

toFront :: Plane -> Vector -> Vector
toFront (plPoint, vPerp) point = if dot point vPerp < dot plPoint vPerp then planeReflect (plPoint, vPerp) point else point

atFront :: Plane -> Vector -> Int       -- 0 is no, 1 is yes, 2 is it's on the plane
atFront (plPoint, vPerp) point = if dot point vPerp > dot plPoint vPerp + 0.000000001       -- anything on the plane is seen as at the front
                                 then 1
                                 else if dot point vPerp < dot plPoint vPerp - 0.000000001
                                      then 0
                                      else 2

translatePlane :: Vector -> Plane -> Plane
translatePlane v (plPoint, vPerp) = (plPoint +. v, vPerp)

rotatePlane :: Quaternion -> Plane -> Plane
rotatePlane q (plPoint, vPerp) = (rotateVect q plPoint, rotateVect q vPerp)

scalePlane :: Double -> Plane -> Plane
scalePlane s (plPoint, vPerp) = (s *. plPoint, s *. vPerp)

trianglePlane :: Triangle -> Plane
trianglePlane (p1,p2,p3) = (p1, cross (50 *. (p2+.((-1)*.p1))) (50 *. (p3+.((-1)*.p1))))

triangleReflect :: Plane -> Triangle -> Triangle
triangleReflect plane (p1,p2,p3) = let r = planeReflect plane in (r p1, r p2, r p3)

translateTriangle :: Vector -> Triangle -> Triangle
translateTriangle offset (p1,p2,p3) = let t = (\v->v+.offset) in (t p1, t p2, t p3)

rotateTri :: Quaternion -> Triangle -> Triangle
rotateTri orient (p1,p2,p3) = let r = rotateVect orient in (r p1, r p2, r p3)

scaleTri :: Double -> Triangle -> Triangle
scaleTri factor (p1,p2,p3) = let s = (\v->factor*.v) in (s p1, s p2, s p3)

triangleVPerp :: Triangle -> Vector
triangleVPerp (p1,p2,p3) = cross (50 *. (p2+.((-1)*.p1))) (50 *. (p3+.((-1)*.p1)))

(+..) :: Vect2d -> Vect2d -> Vect2d
(a,b) +.. (c,d) = (a+c, b+d)

(*..) :: Double -> Vect2d -> Vect2d
c *.. (a,b) = (c*a, c*b)

dot2d :: Vect2d -> Vect2d -> Double
dot2d (a,b) (c,d) = a*c + b*d

project2d :: Vect2d -> Vect2d -> Vect2d
project2d v w = ((dot2d v w)/(dot2d v v)) *.. v

rotateVect2d :: Double -> Vect2d -> Vect2d
rotateVect2d angle (x,y) = ((cos angle) * x - (sin angle) * y, (sin angle) * x + (cos angle) * y)

getAngle :: Vect2d -> Double
getAngle (x,y) = atan2 y x

lineSegmentDist :: Matrix -> Line -> Vect2d -> Double    -- this function cannot assume an orthonormal basis
lineSegmentDist metric (p1,p2) point = let direction = p2 +.. ((-1) *.. p1)
                                           pointDisp = (point +.. ((-1) *.. p1))
                                       in if dot2d direction (matVecMult metric point) <= dot2d direction (matVecMult metric p1)
                                          then (\v -> dot2d v (matVecMult metric v)) (point +.. ((-1) *.. p1))
                                          else if dot2d direction (matVecMult metric point) >= dot2d direction (matVecMult metric p2)
                                          then (\v -> dot2d v (matVecMult metric v)) (point +.. ((-1) *.. p2))
                                          else let projPoint = (\v w -> ((dot2d v (matVecMult metric w)) / (dot2d v (matVecMult metric v))) *.. v) direction pointDisp
                                               in (\v -> dot2d v (matVecMult metric v)) (pointDisp +.. ((-1) *.. projPoint))

triangleDist :: Vector -> Triangle -> Double
triangleDist point (p1,p2,p3) = sqrt( sqPerpDist + sqParrDist )
  where pointDisplacement = point +. ((-1) *. p1)
        e1 = p2 +. ((-1) *. p1)
        e2 = p3 +. ((-1) *. p1)
        metric = ((dot e1 e1, dot e1 e2), (dot e2 e1, dot e2 e2))
        invMetric = matInverse metric
        coPointProjected = (dot e1 pointDisplacement, dot e2 pointDisplacement)
        pointProjected = matVecMult invMetric coPointProjected
        perpVect = project (cross e1 e2) pointDisplacement
        sqPerpDist = (\v -> dot v v) perpVect
        sqParrDist = if fst pointProjected >= 0 && snd pointProjected >= 0 && fst pointProjected + snd pointProjected <= 1
                     then 0
                     else let distL1 = lineSegmentDist metric ((0,1),(1,0)) pointProjected
                              distL2 = lineSegmentDist metric ((0,0),(0,1)) pointProjected
                              distL3 = lineSegmentDist metric ((1,0),(0,0)) pointProjected
                          in min (min distL1 distL2) distL3

triTruncate :: Plane -> Triangle -> [Triangle]
triTruncate plane (p1,p2,p3) = let fronts = [point | point <- [p1,p2,p3], atFront plane point == 1]
                                   backs  = [point | point <- [p1,p2,p3], atFront plane point == 0]
                                   ons    = [point | point <- [p1,p2,p3], atFront plane point == 2]
                               in      if length fronts + length ons == 3 && length ons /= 3  -- one unobscured triangle
                                  then [(p1,p2,p3)]
                                  else if length fronts == 2  -- two triangles make a quad
                                  then let front1 = fronts !! 0
                                           front2 = fronts !! 1
                                           back = backs !! 0
                                           intersect1 = ( ( factorToMeetPlane (translatePlane ((-1) *. front1) plane) (back +. ((-1) *. front1)) ) *. (back +. ((-1) *. front1)) ) +. front1
                                           intersect2 = ( ( factorToMeetPlane (translatePlane ((-1) *. front2) plane) (back +. ((-1) *. front2)) ) *. (back +. ((-1) *. front2)) ) +. front2
                                       in [(front1, intersect1, intersect2), (front2, front1, intersect2)]
                                  else if length fronts == 1 && length ons == 1
                                  then let front = fronts !! 0
                                           back = backs !! 0
                                           on = ons !! 0
                                           intersect = ( ( factorToMeetPlane (translatePlane ((-1) *. front) plane) (back +. ((-1) *. front)) ) *. (back +. ((-1) *. front)) ) +. front
                                       in [(front, on, intersect)]
                                  else if length fronts == 1  -- one triangle, cut off at the plane
                                  then let front = fronts !! 0
                                           back1 = backs !! 0
                                           back2 = backs !! 1
                                           intersect1 = ( ( factorToMeetPlane (translatePlane ((-1) *. front) plane) (back1 +. ((-1) *. front)) ) *. (back1 +. ((-1) *. front)) ) +. front
                                           intersect2 = ( ( factorToMeetPlane (translatePlane ((-1) *. front) plane) (back2 +. ((-1) *. front)) ) *. (back2 +. ((-1) *. front)) ) +. front
                                       in [(front, intersect1, intersect2)]
                                  else []                     -- triangle fully behind plane. Nothing to see at all.


distance :: Vector -> Scene -> Double
distance camPos (Mirror scene mPlane) = distance (toFront mPlane camPos) (ignoreBack mPlane scene)
  
  where ignoreBack plane (Mirror    scene0 plane0) = Mirror (ignoreBack plane0 scene0) plane   -- my reasoning for this is difficult to explain but hopefully will work
        ignoreBack plane (Reflect   scene0 plane0) = Reflect (ignoreBack (planePlaneReflect plane0 plane) scene0) plane0
        ignoreBack plane (Translate scene0 offset) = Translate (ignoreBack (translatePlane ((-1) *. offset) plane) scene0) offset
        ignoreBack plane (Rotate    scene0 orient) = Rotate (ignoreBack (rotatePlane (inv orient) plane) scene0) orient
        ignoreBack plane (ScaleO    scene0 factor) = ScaleO (ignoreBack (scalePlane factor plane) scene0) factor
        ignoreBack plane (Triangles []) = Triangles []
        ignoreBack plane (Triangles (triangle:triangles)) = Triangles ( (triTruncate plane triangle) ++ ((\(Triangles x) -> x) (ignoreBack plane (Triangles triangles))) )


distance camPos (Reflect    scene rPlane) = distance (planeReflect rPlane camPos) scene
distance camPos (Translate  scene offset) = distance (camPos +. ((-1) *. offset)) scene
distance camPos (Rotate     scene orient) = distance (rotateVect (inv orient) camPos) scene
distance camPos (ScaleO     scene factor) = factor * (distance ((1/factor) *. camPos) scene)
distance camPos (Triangles       []     ) = 1/0     --Infinity > all
distance camPos (Triangles (triangle:ts)) = min (triangleDist camPos triangle) (distance camPos (Triangles ts))


nearestTriVPerp :: Vector -> Scene -> Vector
nearestTriVPerp camPos (Mirror scene mPlane) = let rawVPerp = nearestTriVPerp (toFront mPlane camPos) (ignoreBack mPlane scene)
                                               in if atFront mPlane camPos > 0
                                                  then rawVPerp
                                                  else planeReflect mPlane rawVPerp

  where ignoreBack plane (Mirror    scene0 plane0) = Mirror (ignoreBack plane0 scene0) plane
        ignoreBack plane (Reflect   scene0 plane0) = Reflect (ignoreBack (planePlaneReflect plane0 plane) scene0) plane0
        ignoreBack plane (Translate scene0 offset) = Translate (ignoreBack (translatePlane ((-1) *. offset) plane) scene0) offset
        ignoreBack plane (Rotate    scene0 orient) = Rotate (ignoreBack (rotatePlane (inv orient) plane) scene0) orient
        ignoreBack plane (ScaleO    scene0 factor) = ScaleO (ignoreBack (scalePlane factor plane) scene0) factor
        ignoreBack plane (Triangles []) = Triangles []
        ignoreBack plane (Triangles (triangle:triangles)) = Triangles ( (triTruncate plane triangle) ++ ((\(Triangles x) -> x) (ignoreBack plane (Triangles triangles))) )

nearestTriVPerp camPos (Reflect   scene rPlane) = planeReflect ((\(plPoint, vPerp) -> ((0,0,0), vPerp)) rPlane) (nearestTriVPerp (planeReflect rPlane camPos) scene)
nearestTriVPerp camPos (Translate scene offset) = nearestTriVPerp (camPos +. ((-1) *. offset)) scene
nearestTriVPerp camPos (Rotate    scene orient) = rotateVect orient (nearestTriVPerp (rotateVect (inv orient) camPos) scene)
nearestTriVPerp camPos (ScaleO    scene factor) = factor *. (nearestTriVPerp ((1/factor) *. camPos) scene)
nearestTriVPerp camPos (Triangles triangles   ) = let closest = fst (nearestTriangle triangles)
                                                      rawVPerp = triangleVPerp closest
                                                      signedVPerp = if dot rawVPerp (camPos +. ((-1) *. ((\(a,b,c) -> a) closest))) < 0
                                                                    then (-1) *. rawVPerp
                                                                    else rawVPerp
                                                  in (\v -> (1 / (sqrt (dot v v)) ) *. v) signedVPerp
  where nearestTriangle (triangle:[]) = (triangle, triangleDist camPos triangle)
        nearestTriangle (triangle:ts) = let prevOutput = nearestTriangle ts
                                        in if triangleDist camPos triangle < snd prevOutput
                                           then (triangle, triangleDist camPos triangle)
                                           else prevOutput
        nearestTriangle [] = error "Failed to find a triangle."


--compressScene :: Scene -> CompressedScene
--compressScene (Mirror scene plane) = ignoreBack ( compressScene scene ) plane
--  where ignoreBack ([]    , mirrorPlanes) plane = ([], (plane:mirrorPlanes) )
--        ignoreBack ((t:ts), mirrorPlanes) plane = let (prevTs, prevMirrorPlanes) = ignoreBack (ts, mirrorPlanes) plane
--                                                  in ((triTruncate plane t) ++ prevTs, prevMirrorPlanes)

compressScene (Mirror    scene mirrorPlane) = (\(a,b)->(a,mirrorPlane:b)) $ compressScene scene

compressScene (Reflect   scene rPlane) = reflectCompSc (compressScene scene) rPlane
  where reflectCompSc ([]    , []            ) rPlane = ([]                               , []                                         )
        reflectCompSc ([]    , (plane:planes)) rPlane = ([]                               , (planePlaneReflect rPlane plane):prevPlanes)
          where (_     , prevPlanes) = reflectCompSc ([], planes) rPlane
        reflectCompSc ((t:ts), planes        ) rPlane = ((triangleReflect rPlane t):prevTs, prevPlanes                                 )
          where (prevTs, prevPlanes) = reflectCompSc (ts, planes) rPlane

compressScene (Translate scene offset) = translateCompSc (compressScene scene) offset
  where translateCompSc ([]    , []            ) offset = ([]                                 , []                                      )
        translateCompSc ([]    , (plane:planes)) offset = ([]                                 , (translatePlane offset plane):prevPlanes)
          where (_     , prevPlanes) = translateCompSc ([], planes) offset
        translateCompSc ((t:ts), planes        ) offset = ((translateTriangle offset t):prevTs, prevPlanes                              )
          where (prevTs, prevPlanes) = translateCompSc (ts, planes) offset

compressScene (Rotate    scene orient) = rotateCompSc (compressScene scene) orient
  where rotateCompSc ([]    , []            ) orient = ([]                         , []                                   )
        rotateCompSc ([]    , (plane:planes)) orient = ([]                         , (rotatePlane orient plane):prevPlanes)
          where (_     , prevPlanes) = rotateCompSc ([], planes) orient
        rotateCompSc ((t:ts), planes        ) orient = ((rotateTri orient t):prevTs, prevPlanes                           )
          where (prevTs, prevPlanes) = rotateCompSc (ts, planes) orient

compressScene (ScaleO    scene factor) = scaleCompSc (compressScene scene) factor
  where scaleCompSc ([]    , []            ) factor = ([]                        , []                                  )
        scaleCompSc ([]    , (plane:planes)) factor = ([]                        , (scalePlane factor plane):prevPlanes)
          where (_     , prevPlanes) = scaleCompSc ([], planes) factor
        scaleCompSc ((t:ts), planes        ) factor = ((scaleTri factor t):prevTs, prevPlanes                          )
          where (prevTs, prevPlanes) = scaleCompSc (ts, planes) factor

compressScene (Triangles ts) = (ts, [])



compScDist :: Vector -> CompressedScene -> Double
compScDist point (tris    , (plane:planes)) = let newPoint = toFront plane point
                                              in  compScDist newPoint (tris, planes)
compScDist point (tri:tris, []            ) = let prevMinDist = compScDist point (tris, [])
                                              in  min (triangleDist point tri) prevMinDist
compScDist point ([]      , []            ) = 1/0



compNearestTriVPerp :: Vector -> CompressedScene -> Vector
compNearestTriVPerp vect (tris, (plane:planes)) = if atFront plane vect /= 0
                                                  then compNearestTriVPerp vect (tris, planes)
                                                  else planeReflect ((\(plPoint, vPerp)->((0,0,0), vPerp)) plane) $ compNearestTriVPerp (planeReflect plane vect) (tris, planes)
compNearestTriVPerp vect (tris, []            ) = let closest = fst (nearestTriangle tris)
                                                      rawVPerp = triangleVPerp closest
                                                      signedVPerp = if dot rawVPerp (vect +. ((-1) *. ((\(a,b,c) -> a) closest))) < 0
                                                                    then (-1) *. rawVPerp
                                                                    else rawVPerp
                                                  in (\v -> (1 / (sqrt (dot v v)) ) *. v) signedVPerp
  where nearestTriangle (triangle:[]) = (triangle, triangleDist vect triangle)
        nearestTriangle (triangle:ts) = let prevOutput = nearestTriangle ts
                                        in if triangleDist vect triangle < snd prevOutput
                                           then (triangle, triangleDist vect triangle)
                                           else prevOutput
        nearestTriangle [] = error "Failed to find a triangle."


marchRayToScene :: CompressedScene -> Double -> Double -> (Vector, Vector) -> (Bool, Integer, Double, Vector)      -- output is (hitScene, stepsToHit, closestPass, rayPos). Not all will be meaningful. Also, ray is a unit ray.
marchRayToScene scene hitThreshold allowedTravel ray = let camDist = 0.99 * (compScDist (fst ray) scene)      --the 0.99 is to make sure it never goes past. This can be an issue with an exact distance estimator.
                                                       in if camDist > allowedTravel
                                                          then (False, 0, camDist, (0,0,0))
                                                          else if camDist < hitThreshold
                                                          then (True, 0, 0.0, fst ray)
                                                          else let (hts, sth, clp, r) = marchRayToScene scene hitThreshold (allowedTravel - camDist) ((fst ray +. (camDist *. (snd ray))), snd ray)
                                                               in (hts, sth+1, min clp camDist, r)


marchRayToLight :: CompressedScene -> Double -> Double -> Vector -> Vector -> Double -> Double        --Output is min angle
marchRayToLight scene hitThreshold allowedTravel lightDir rayPos distToNow = let rayDist = 0.99 * (compScDist rayPos scene)
                                                                                 newDistToNow = distToNow + rayDist
                                                                                 currentAngle = angleFromSphere rayDist distToNow
                                                                             in if distToNow > allowedTravel
                                                                                then currentAngle
                                                                                else if rayDist < hitThreshold
                                                                                then 0
                                                                                else let prevAngle = marchRayToLight scene hitThreshold allowedTravel lightDir (rayPos +. (rayDist *. lightDir)) newDistToNow
                                                                                     in  min prevAngle currentAngle

angleFromSphere :: Double -> Double -> Double
angleFromSphere r d = if r > d                                             --it can be assumed that r is positive
                      then error "turns out this case is real"             --Not the actual angle. This input is mostly invalid but I'm suspicious that "almost equal" shenanigans might make this necessary.
                      else asin (r/d)


constructRay :: Int -> Int -> Double -> (Vector, Vector)
--constructRay x y = ( ((fromIntegral x) * 0.01 - 1.5, 2 - (fromIntegral y) * 0.01, -2), (0.0, 0.0, 1.0) )                 --This creates parallel projection
constructRay x y fov = let phi = atan (sqrt ((fov * (fromIntegral x - 960) )^2 + (fov * (fromIntegral y - 540) )^2))     --This creates perspective projection
                           theta = atan2 (fov * (fromIntegral y - 540) ) (fov * (fromIntegral x - 960) )
                       in ((0,0,-5), (cos theta * sin phi, -sin theta * sin phi, cos phi))


getPixelColour :: CompressedScene -> Vector -> Double -> Int -> Int -> PixelRGB8
getPixelColour scene lightDir fov x y = let (hitScene, stepsTaken, closestPass, hitPos) = marchRayToScene scene 0.01 10.0 (constructRay x y fov)
                                    in if hitScene
                                       then let brightness = truncate (abs (sqrt ( ambientLight +  sourceLight)))
                                                ambientLight = 0.8 ^ stepsTaken * 9000
                                                facingSourceFactor = dot (compNearestTriVPerp hitPos scene) lightDir
                                                angleFromLightMarch = marchRayToLight scene ((compScDist hitPos scene) * 0.99) 10.0 lightDir hitPos (compScDist hitPos scene)
                                                occludedSourceFactor = 1--min (angleFromLightMarch/0.1) 1
                                                sourceLight  = 45000 * facingSourceFactor * occludedSourceFactor
                                            in PixelRGB8 brightness brightness brightness
                                       else PixelRGB8 0 200 255


sierpTetra :: Integer -> Scene
sierpTetra n | n == 0 = let p1 = ( 1, 1, 1)
                            p2 = ( 1,-1,-1)
                            p3 = (-1, 1,-1)
                            p4 = (-1,-1, 1)
                        in Triangles [(p1,p2,p3), (p1,p2,p4), (p1,p3,p4), (p2,p3,p4)]
             | n >  0 = Mirror (
                          Mirror (
                            Mirror (
                              Translate ( ScaleO ( Translate ( sierpTetra (n-1) ) (-1,-1,-1) ) 0.5 ) (1,1,1)
                              )
                            ((0,0,0),(0,1,1))
                            )
                          ((0,0,0),(1,0,1))
                          )
                        ((0,0,0),(1,1,0))

compSierp :: Integer -> CompressedScene
compSierp n = compressScene ( Rotate (sierpTetra n) $ qMult (constructRotation (0,1,0) t) (constructRotation (cross ((\v->(1/(sqrt $ dot v v)) *. v) (1,1,1)) (0,1,0)) (dot ((\v->(1/(sqrt $ dot v v)) *. v) (1,1,1)) (0,1,0)) ) )
  where t = 8*pi/12

sierpImage = generateImage (getPixelColour ( compSierp 7 ) ((\v->(1/(sqrt $ dot v v)) *. v) (-0.3, 1, -1)) 0.001) 1920 1080

triFrac :: Integer -> Scene
triFrac n | n == 0 = Triangles [((0,0,1),(-(sqrt 3)/2,0,-1/2),(sqrt(3)/2,0,-1/2))]
          | n >  0 = Mirror (
                     Mirror (
                     Mirror (
                            Translate ( ScaleO ( triFrac (n-1) ) 0.5 ) (0,0,0.5)
                            ) ((0,0,1/4),(0,-sqrt(3),1))
                            ) ((0,0,0  ),( 1,0,sqrt(3)))
                            ) ((0,0,0  ),(-1,0,sqrt(3)))

compTriFrac :: Integer -> CompressedScene
compTriFrac n = compressScene ( Rotate (ScaleO(triFrac n) 1.8) $ qMult (constructRotation (0,1,0) t) (constructRotation (-1,0,0) s) )
  where t = -1*pi/12
        s = 1*pi/6

triFracImage = generateImage (getPixelColour ( compTriFrac 7 ) ((\v->(1/(sqrt $ dot v v)) *. v) (-1, 1, -0.8)) 0.001) 1920 1080

sponge :: Integer -> Scene
sponge n | n == 0 = Mirror (
                    Mirror (
                    Mirror (
                    Mirror (
                    Mirror (
                           Triangles [((1.5,1.5,1.5),(-1.5,1.5,1.5),(1.5,-1.5,1.5))]
                           ) ((0,0,0),(0,-1,1))
                           ) ((0,0,0),(-1,0,1))
                           ) ((0,0,0),(1,1,0))
                           ) ((0,0,0),(1,0,1))
                           ) ((0,0,0),(0,1,1))
         | n > 0  = Mirror (
                    Mirror (
                    Mirror (
                    Mirror (
                    Mirror (
                    Mirror (
                           Translate (ScaleO (sponge (n-1)) (1/3)) (1,1,1)
                           ) ((0.5,0,0),(1,0,0))
                           ) ((0,0,0),(1,0,0))
                           ) ((0,0,0),(-1,1,0))
                           ) ((0,0,0),(1,1,0))
                           ) ((0,0,0),(-1,0,1))
                           ) ((0,0,0),(1,0,1))

compSponge :: Integer -> CompressedScene
compSponge n = compressScene (Rotate (sponge n) $ constructRotation (0,1,0) t)
  where t = 3*pi/8

spongeImage = generateImage (getPixelColour ( compSponge 7 ) ((\v->(1/(sqrt $ dot v v)) *. v) (-1, 3, -0.8)) 0.001) 1920 1080

writeImage name image = writePng name image

main = writeImage "sponge.png" spongeImage
