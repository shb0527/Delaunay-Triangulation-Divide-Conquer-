#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <math.h>
#include <limits>
#include <tuple>
#include <type_traits>
#include <SFML/Graphics.hpp>
namespace krs {



    template <typename T>
    struct Vector2
    {
        using Type = T;
        Vector2() = default;
        Vector2(const Vector2<T>& v) = default;
        Vector2(const T vx, const T vy);
        template<typename U>
        friend std::ostream& operator <<(std::ostream& str, const Vector2<U>& v);
        T x;
        T y;

        T norm2() const;
        T dist2(const Vector2<T>& v) const;
        T dist(const Vector2<T>& v) const;
        Vector2& operator=(const Vector2<T>&) = default;
        Vector2& operator=(Vector2&&) = default;
        bool operator == (const Vector2<T>& v) const;

    };

    template <typename T>
    Vector2<T>::Vector2(const T vx, const T vy) :
        x(vx), y(vy)
    {}

    template<typename T>
    bool
        Vector2<T>::operator ==(const Vector2<T>& v) const
    {
        return (this->x == v.x) && (this->y == v.y);
    }
    template<typename U>
    std::ostream&
        operator <<(std::ostream& str, const Vector2<U>& v)
    {
        return str << "Point x: " << v.x << " y: " << v.y;
    }

    template<typename T>
    T
        Vector2<T>::norm2() const
    {
        return x * x + y * y;
    };

    template<typename T>
    T
        Vector2<T>::dist2(const Vector2<T>& v) const
    {
        const T dx = x - v.x;
        const T dy = y - v.y;
        return dx * dx + dy * dy;
    };
    template<>
    double
        Vector2<double>::dist(const Vector2<double>& v) const
    {
        return hypot(x - v.x, y - v.y);
    };


    template<class T>
    typename std::enable_if<std::is_same <T, double>::value, bool> ::type
        almost_equal(T x, T y, int ulp = 2)
    {
        return fabsf(x - y) <= std::numeric_limits<double>::epsilon() * fabsf(x + y) * static_cast<double>(ulp)
            || fabsf(x - y) < std::numeric_limits<double>::min();

    }
    template<typename T>
    bool almost_equal(const Vector2<T>& v1, const Vector2<T>& v2)
    {
        return almost_equal(v1.x, v2.x) && almost_equal(v1.y, v2.y);
    }




    template<typename T>
    bool almost_equal2(std::pair<T, T>& p1, std::pair<T, T>& p2)
    {
        return almost_equal(p1.first, p2.first) && almost_equal(p1.second, p2.second);
    }

    template<typename T>
    bool containsonePoint(std::pair<std::pair<T, T>, std::pair<T, T> > s, std::pair<T, T> p1, std::pair<T, T> p2)
    {

        return (almost_equal2(s.first, p1) || almost_equal2(s.second, p2)) || ((almost_equal2(s.first, p2) || almost_equal2(s.second, p1)));

    }
    template<typename T>
    bool containstwoPoints(std::pair<std::pair<T, T>, std::pair<T, T> > s, std::pair<T, T> p1, std::pair<T, T> p2)
        {

            return (almost_equal2(s.first, p1) || almost_equal2(s.second, p2)) && ((almost_equal2(s.first, p2) || almost_equal2(s.second, p1)));
        }
    template<typename T>
    struct Triangle {
        using Type = T;
        using VertexType = Vector2<Type>;
   

        Triangle() = default;
        Triangle(const Triangle&) = default;
        Triangle(Triangle&&) = default;
        Triangle(const VertexType& v1, const VertexType& v2, const VertexType& v3);

        template<typename U>
        friend std::ostream& operator <<(std::ostream& str, const Triangle<U>& t);

        bool circumCircleContains(const VertexType& v) const;



        const VertexType* a;
        const VertexType* b;
        const VertexType* c;

    };



    template<typename T>
    Triangle<T>::Triangle(const VertexType& v1, const VertexType& v2, const VertexType& v3) :
        a(&v1), b(&v2), c(&v3)
    {}
    template<typename U>
    std::ostream&
        operator <<(std::ostream& str, const Triangle<U>& t)
    {
        return str << "Triangle:" << "\n\t" <<
            *t.a << "\n\t" <<
            *t.b << "\n\t" <<
            *t.c << '\n';
    }
    template<typename T>
    bool
        Triangle<T>::circumCircleContains(const Vector2<T>& v) const
    {
        const T ab = a->norm2();
        const T cd = b->norm2();
        const T ef = c->norm2();

        const T ax = a->x;
        const T ay = a->y;
        const T bx = b->x;
        const T by = b->y;
        const T cx = c->x;
        const T cy = c->y;

        const T circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
        const T circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));

        const VertexType circum(circum_x / 2, circum_y / 2);
        const T circum_radius = a->dist2(circum);
        const T dist = v.dist2(circum);


        return dist <= circum_radius;
    };

    template <typename T>
    class delaunay {
        using VertexType = Vector2<T>;

        using Type = T;
        using TriangleType = Triangle<T>;

        //std::vector<Edge<double>> _edges;
        std::vector<Vector2<T>> _vertices;
        std::vector <Triangle<T>> _triangles;
        std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >> _edges;

    public:

        delaunay() = default;
        delaunay(const delaunay&) = default;
        delaunay(delaunay&&) = default;


        const std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >> potentials(std::vector<VertexType>& vertices1, std::vector<VertexType>& vertices2
        ,  std::pair<T, T> new_pair1, std::pair<T, T> new_pair2, bool isStart, int depth);
        void partition(std::vector<VertexType>& vertices, size_t begin, size_t end);
        const std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >> getEdges();
        const std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >> createEdges(std::vector<VertexType>& vertices);
        const std::vector<std::pair<std::pair<T, T>, std::pair<T, T>>> slicingVector(std::vector<VertexType>& vertices, size_t b, size_t e);
       
        bool onsegment(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const;
        const int orientation(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const;
        bool intersect(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> q1, std::pair<T, T> q2) const;
        bool delaunay_intersect(std::pair<T, T> np1, std::pair<T, T> np2) const;
        delaunay& operator=(const delaunay&) = default;
        delaunay& operator=(delaunay&&) = default;

    };



    template<typename T>
    const std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >>
        delaunay<T>::getEdges()
    {
        return _edges;
    }


    

   template<typename T>
    const std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >>
        delaunay<T>::slicingVector(std::vector<Vector2<T>>& vertices, size_t b, size_t e)
    {
        std::pair<T, T> np1 = { 0,0 };
        std::pair<T, T> np2 = { 0,0 };
        std::vector<std::vector<Vector2<T>>> vlist;
        for (auto i = 0; i < vertices.size(); i += 3) {
            auto start = vertices.begin() + i;
            auto end = vertices.begin() + i + 3;

            std::vector<Vector2<T>> v(start, end);
            vlist.push_back(v);
            createEdges(v);
        }   
       
        
        for (auto i = 0; i < vlist.size()-1; i++)
        {
            potentials(vlist[i], vlist[i + 1], np1, np2, true, 0);
        }


      std::vector<std::vector<Vector2<T>>> vlist2;
        for (auto i = 0; i < vertices.size(); i += 4) {
             auto start = vertices.begin() + i;
             auto end = vertices.begin() + i + 4;

             std::vector<Vector2<T>> v(start, end);
             vlist2.push_back(v);
         }
        for (auto i = 0; i < vlist2.size() - 1; i++)
        {
            potentials(vlist2[i], vlist2[i + 1], np1, np2, true, 0);
        }

         std::vector<std::vector<Vector2<T>>> vlist3;
         for (auto i = 0; i < vertices.size(); i += 6) {
             auto start = vertices.begin() + i;
             auto end = vertices.begin() + i + 6;

             std::vector<Vector2<T>> v(start, end);
             vlist3.push_back(v);
         }
        
         potentials(vlist3[0], vlist3[1], np1, np2, true, 0);
         potentials(vlist3[1], vlist3[2], np1, np2, true, 0);
         potentials(vlist3[2], vlist3[3], np1, np2, true, 0);

         std::vector<std::vector<Vector2<T>>> vlist4;
         for (auto i = 0; i < vertices.size(); i += 8) {
             auto start = vertices.begin() + i;
             auto end = vertices.begin() + i + 8;

             std::vector<Vector2<T>> v(start, end);
             vlist4.push_back(v);
         }
         
         potentials(vlist4[0], vlist4[1], np1, np2, true, 0);
         potentials(vlist4[1], vlist4[2], np1, np2, true, 0);

         std::vector<std::vector<Vector2<T>>> vlist5;
         for (auto i = 0; i < vertices.size(); i += 12) {
             auto start = vertices.begin() + i;
             auto end = vertices.begin() + i + 12;

             std::vector<Vector2<T>> v(start, end);
             vlist5.push_back(v);
         }
   
         potentials(vlist5[0], vlist5[1], np1, np2, true, 0);
           
         return _edges;
     
    }
    template <typename T>
    const std::vector<std::pair<std::pair<T,T>, std::pair<T,T> >> 
        delaunay<T>::createEdges(std::vector<Vector2<T>>& vertices) {

        
        size_t s1 = vertices.size() - 1;
        
        if (s1 == 1) {

            T firstx = vertices[0].x;
            T firsty = vertices[0].y;

       
            std::pair<T, T> pair1 = { firstx, firsty };
            std::pair<T, T> pair2 = { vertices[s1].x , vertices[s1].y };
           
            _edges.push_back({ pair1, pair2 });


        }
        else if (s1 == 2) {


            T firstx = vertices[0].x;
            T firsty = vertices[0].y;

          
            const std::pair<double, double> pair1 = { firstx, firsty };
            const std::pair<double, double> pair2 = { vertices[s1-1].x , vertices[s1-1].y };
            const std::pair<double, double> pair3 = { vertices[s1].x , vertices[s1].y };
            
            _edges.push_back({ pair1, pair2 });
            _edges.push_back({ pair2, pair3 });
            _edges.push_back({ pair3, pair1 });

        };
        return _edges;
    }
    template<typename T>
    bool delaunay<T>::onsegment(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const
    {
        if (p2.first <= std::max(p1.first, p3.first) && p2.first >= std::min(p1.first, p3.first)
            && p2.second <= std::max(p1.second, p3.second) && p2.second >= std::min(p1.second, p3.second) ){
            return true;
        }
        return false;
    }

    template<typename T>
    const int delaunay<T>::orientation(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> p3) const {

        T val = (p2.second - p1.second) * (p3.first - p2.first) -
            (p2.first - p1.first) * (p3.second - p2.second);
        if (val == 0) {
            return 0;
        }
        return (val > 0) ? 1 : 2;
    }

    template<typename T>
    bool delaunay<T>::intersect(std::pair<T, T> p1, std::pair<T, T> p2, std::pair<T, T> q1, std::pair<T, T> q2) const
    {
        T o1 = orientation(p1, q1, p2);
        T o2 = orientation(p1, q1, q2);
        T o3 = orientation(p2, q2, p1);
        T o4 = orientation(p2, q2, q1);

        if (o1 != o2 && o3 != o4)
            return true;

        if (o1 == 0 && onsegment(p1, p2, q1)) return true;
        if (o2 == 0 && onsegment(p1, q2, q1)) return true;
        if (o3 == 0 && onsegment(p2, p1, q2)) return true;
        if (o4 == 0 && onsegment(p2, q1, q2)) return true;

        return false; 

    }

    template<typename T>
    bool delaunay<T>::delaunay_intersect(std::pair<T, T> np1, std::pair<T, T> np2) const {

        int c = 0;
        for (auto& e : _edges) {
            if (intersect(np1, e.first, np2, e.second) && !containsonePoint(e, np1, np2)){ 
                return false;
           }
 
        }return true;
    }


    template <typename T>
    const std::vector<std::pair<std::pair<T, T>, std::pair<T, T> >>
        delaunay<T>::potentials(std::vector<Vector2<T>>& vertices1, std::vector<Vector2<T>>& vertices2,
            std::pair<T, T> np1, std::pair<T, T> np2, bool isStart, int depth) {

      
        if (isStart) {

            std::sort(vertices1.begin(), vertices1.end(), [](Vector2<T> a, Vector2<T> b) { return a.y > b.y; });
            std::sort(vertices2.begin(), vertices2.end(), [](Vector2<T> a, Vector2<T> b) { return a.y > b.y; });

            const VertexType p1(vertices1[0].x, vertices1[0].y);
            const VertexType p2(vertices2[0].x, vertices2[0].y);

            const std::pair<T, T> np1 = { vertices1[0].x, vertices1[0].y };
            const std::pair<T, T> np2 = { vertices2[0].x, vertices2[0].y };
            std::pair<std::pair<T, T>, std::pair<T, T> > enp = { np1, np2 };
            if (delaunay_intersect(np1, np2)) {
                _edges.push_back({ np1, np2 });
            }
            
        
          


            ///******************  RR edge  *************************/

            T start = 0;
            T angle = 0;
            std::vector<Vector2<double>>::iterator it;
            std::vector<Vector2<double>> Right_klist;
            bool isright = false;
            int rcount = 0;
            for (auto& k : vertices2)
            {
                start++;

                std::pair<T, T> p3 = { k.x, k.y };

                double dot = k.x * (p2.x) + k.y * p2.y;
                double ma = sqrt(k.norm2());
                double mb = sqrt(p2.norm2());
                angle = acos(dot / (ma * mb));
                angle = angle * 180 / 3.1415;

                if (angle > 0 && angle < 180)
                {

                    const TriangleType t = TriangleType{ p1, p2, k };
                    int c = 0;
                    for (it = vertices2.begin() + start; it != vertices2.end(); it++) {

                        if (t.circumCircleContains(*it)) {

                            _edges.erase(std::remove_if(begin(_edges), end(_edges), [np2, p3]
                            (std::pair<std::pair<T, T>, std::pair<T, T> >& j) {
                                    return containstwoPoints(j, np2, p3); }), end(_edges));
                            break;

                        }

                        else { c++; }
                    }
                    if (c == (vertices2.size() - start)) {
                        // potential is submitted
                        isright = true;
                        Right_klist.push_back(k);
                        rcount++;
                    }

                }
            }

            ///******************  LL edge  *************************/
            start = 0;
            bool isleft = false;
            int lcount = 0;
            std::vector<Vector2<double>> Left_klist;
            for (auto& k : vertices1)
            {
                start++;

                std::pair<T, T> p3 = { k.x, k.y };

                double dot = k.x * (p1.x) + k.y * p1.y;
                double ma = sqrt(k.norm2());
                double mb = sqrt(p1.norm2());
                angle = acos(dot / (ma * mb));
                angle = angle * 180 / 3.1415;

                if (angle > 0 && angle < 180)
                {

                    const TriangleType t = TriangleType{ p1, p2, k };
                    int c = 0;
                    for (it = vertices1.begin() + start; it != vertices1.end(); it++) {

                        if (t.circumCircleContains(*it)) {
                            _edges.erase(std::remove_if(begin(_edges), end(_edges), [np1, p3]
                            (std::pair<std::pair<T, T>, std::pair<T, T> >& j) {
                                    return containstwoPoints(j, np1, p3); }), end(_edges));
                                    break; }
                        else { c++; }
                    }
                    if (c == (vertices1.size() - start)) {
                        // potential is submitted
                        isleft = true;
                        Left_klist.push_back(k);
                        lcount++;

                    }

                }

            }
           
            if (depth != vertices1.size()) {
                if (isleft == true && isright == false) {
                    std::sort(Left_klist.begin(), Left_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x > b.x; });

                    for (auto& a : Left_klist) {
                        const std::pair<T, T> new_pair1 = { a.x, a.y };
                        const std::pair<T, T> new_pair2 = { vertices2[0].x, vertices2[0].y };
                        if (delaunay_intersect(new_pair1, new_pair2)) {
                            _edges.push_back({ new_pair1, new_pair2 });
                            potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                        }
                    } return _edges;
                }
                else if (isleft == false && isright == true) {
                    std::sort(Right_klist.begin(), Right_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x < b.x; });

                    for (auto& a : Right_klist) {
                        const std::pair<T, T> new_pair1 = { vertices1[0].x, vertices1[0].y };
                        const std::pair<T, T> new_pair2 = { a.x, a.y };
                        if (delaunay_intersect(new_pair1, new_pair2)) {
                            _edges.push_back({ new_pair1, new_pair2 });
                            potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                        }
                    }
                }
                else if (isleft == true && isright == true) {
                    std::sort(Left_klist.begin(), Left_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x > b.x; });
                    std::sort(Right_klist.begin(), Right_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x < b.x; });
                    const TriangleType t = TriangleType{ p1, p2, Left_klist[0] };
                    if (!t.circumCircleContains(Right_klist[0])) {
                        std::sort(Right_klist.begin(), Right_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x < b.x; });
                        for (auto& a : Right_klist) {
                            const std::pair<T, T> new_pair1 = { vertices1[0].x, vertices1[0].y };
                            const std::pair<T, T> new_pair2 = { a.x, a.y };
                            if (delaunay_intersect(new_pair1, new_pair2)) {
                                _edges.push_back({ new_pair1, new_pair2 });
                                potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                            }
                        }
                    }
                    else {

                        std::sort(Left_klist.begin(), Left_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x > b.x; });
                        std::sort(Right_klist.begin(), Right_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x < b.x; });
                        const TriangleType t = TriangleType{ p1, p2, Right_klist[0] };
                        if (!t.circumCircleContains(Left_klist[0])) {

                            for (auto& a : Left_klist) {
                                const std::pair<T, T> new_pair1 = { a.x, a.y };
                                const std::pair<T, T> new_pair2 = { vertices2[0].x, vertices2[0].y };
                                if (delaunay_intersect(new_pair1, new_pair2)) {
                                    _edges.push_back({ new_pair1, new_pair2 });
                                    potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                                }
                            }
                        } return _edges;
                    }
                }
            }
            else { return _edges; }


            
        }


            else {


                const VertexType p1 = VertexType{ np1.first, np1.second };
                const VertexType p2 = VertexType{ np2.first, np2.second };
                //RR Edge 

                T start = 0;
                T angle = 0;

                std::vector<Vector2<double>>::iterator it;


                ///******************  RR edge  *************************/

                std::vector<Vector2<double>> Right_klist;
                bool isright = false;
                int rcount = 0;
                for (auto& k : vertices2)
                {
                    start++;

                    std::pair<T, T> p3 = { k.x, k.y };

                    double dot = k.x * (p2.x) + k.y * p2.y;
                    double ma = sqrt(k.norm2());
                    double mb = sqrt(p2.norm2());
                    angle = acos(dot / (ma * mb));
                    angle = angle * 180 / 3.1415;

                    if (angle > 0 && angle < 180)
                    {

                        const TriangleType t = TriangleType{ p1, p2, k };
                        int c = 0;
                        for (it = vertices2.begin() + start; it != vertices2.end(); it++) {

                            if (t.circumCircleContains(*it)) {
                                ;
                                _edges.erase(std::remove_if(begin(_edges), end(_edges), [np2, p3]
                                (std::pair<std::pair<T, T>, std::pair<T, T> >& j) {
                                        return containstwoPoints(j, np2, p3); }), end(_edges));
                                break;

                            }

                            else { c++; }
                        }
                        if (c == (vertices2.size() - start)) {
                            // potential is submitted
                            isright = true;
                            Right_klist.push_back(k);
                            rcount++;
                        }

                    }
                }

                ///******************  LL edge  *************************/
                start = 0;
                bool isleft = false;
                int lcount = 0;
                std::vector<Vector2<double>> Left_klist;
                for (auto& k : vertices1)
                {
                    start++;

                    std::pair<T, T> p3 = { k.x, k.y };

                    double dot = k.x * (p1.x) + k.y * p1.y;
                    double ma = sqrt(k.norm2());
                    double mb = sqrt(p1.norm2());
                    angle = acos(dot / (ma * mb));
                    angle = angle * 180 / 3.1415;

                    if (angle > 0 && angle < 180)
                    {

                        const TriangleType t = TriangleType{ p1, p2, k };
                        int c = 0;
                        for (it = vertices1.begin() + start; it != vertices1.end(); it++) {

                            if (t.circumCircleContains(*it)) {
                              _edges.erase(std::remove_if(begin(_edges), end(_edges), [np1, p3]
                                (std::pair<std::pair<T, T>, std::pair<T, T> >& j) {
                                        return containstwoPoints(j, np1, p3); }), end(_edges));
                                break;

                            }
                            else { c++; }
                        }
                        if (c == (vertices1.size() - start)) {
                            // potential is submitted
                            isleft = true;
                            Left_klist.push_back(k);
                            lcount++;
                        }

                    }
                }


       

               if (depth != vertices1.size()) {
                    if (isleft == true && isright == false) {
                        std::sort(Left_klist.begin(), Left_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x > b.x; });
                        for (auto& a : Left_klist) {

                            const std::pair<T, T> new_pair1 = { a.x, a.y };
                            const std::pair<T, T> new_pair2 = { np2.first, np2.second };

                            if (delaunay_intersect(new_pair1, new_pair2)) { 
                                _edges.push_back({ new_pair1, new_pair2 });
                                potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                            }
                        }
                        return _edges;
                       
                    }
                    else if (isleft == false && isright == true) {
                        std::sort(Right_klist.begin(), Right_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x < b.x; });

                        for (auto &a: Right_klist) {
                            const std::pair<T, T> new_pair1 = { np1.first, np1.second };
                            const std::pair<T, T> new_pair2 = { a.x, a.y };
                            if (delaunay_intersect(new_pair1, new_pair2)) {
                                _edges.push_back({ new_pair1, new_pair2 });
                                potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                            }

                        } return _edges;
                       
                    }
                    else if (isleft == true && isright == true) {

                        std::sort(Left_klist.begin(), Left_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x > b.x; });
                        std::sort(Right_klist.begin(), Right_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x < b.x; });
                        const TriangleType t = TriangleType{ p1, p2, Left_klist[0] };

                        if (!t.circumCircleContains(Right_klist[0])) {

                            for (auto& a : Left_klist) {
                                const std::pair<T, T> new_pair1 = { a.x, a.y };
                                const std::pair<T, T> new_pair2 = { np2.first, np2.second };
                                if (delaunay_intersect(new_pair1, new_pair2)) {
                                    _edges.push_back({ new_pair1, new_pair2 });
                                    potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                                }
                            } 
                        }
                        else {

                            std::sort(Left_klist.begin(), Left_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x > b.x; });
                            std::sort(Right_klist.begin(), Right_klist.end(), [](Vector2<T> a, Vector2<T> b) { return a.x < b.x; });
                            const TriangleType t = TriangleType{ p1, p2, Right_klist[0] };
                            if (!t.circumCircleContains(Left_klist[0])) {

                                for (auto& a : Right_klist) {
                                    const std::pair<T, T> new_pair1 = { np1.first, np1.second };
                                    const std::pair<T, T> new_pair2 = { a.x, a.y };
                                    if (delaunay_intersect(new_pair1, new_pair2)) {
                                        _edges.push_back({ new_pair1, new_pair2 });
                                        potentials(vertices1, vertices2, new_pair1, new_pair2, false, depth + 1);
                                    }

                                }
                            } return _edges;
                        }
                    }
                    else { return _edges; }
                }

            }


            return _edges;
        };

    
}




int main() {


    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<> dist_w(0, 800);
    std::uniform_real_distribution<> dist_h(0, 600);

    std::vector<krs::Vector2<double>> points;

   for (int i = 0; i < 24; ++i) {
        points.push_back(krs::Vector2<double>{dist_w(gen), dist_h(gen)});
    }
    std::sort(points.begin(), points.end(), [](krs::Vector2<double> a, krs::Vector2<double> b) { return a.x < b.x; });
    for (auto& a : points) {
        std::cout << "[" << a.x << ", " << a.y << "] ";
    }



    size_t size = points.size() - 1;
    size_t mid = points.size() / 2;

    krs::delaunay<double> delaunay;


    std::pair<std::pair<double, double>, std::pair<double, double> > pair1;

    std::pair<double, double> pair3 = { 2.3, 4.5 };
    pair1 = { { 2.3, 4.5 } , { 3.2, 3.12} };


    //delaunay.partition(points, 0, size);
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double> >> edges = delaunay.slicingVector(points, 0, size);
    //std::vector<std::pair<std::pair<double, double>, std::pair<double, double> >> edges = delaunay.getEdges();
    std::cout << std::endl;


  /*  for (auto& e : edges) {
        std::cout << "(" << e.first.first << ", " << e.first.second << ")" << " (" << e.second.first << ", " << e.second.second << ")" << std::endl;
    }
    std::cout << edges.size() << std::endl;*/

    std::vector<std::pair<std::pair<double, double>, std::pair<double, double> >> medges;
    
       //double p1 = dist_w(gen);
       // double p2 = dist_w(gen);
       // double p3 = dist_h(gen);
       // double p4 = dist_h(gen);
       // double p5 = dist_w(gen);
       // double p6 = dist_w(gen);
       // double p7 = dist_h(gen);
       // double p8 = dist_h(gen);

       // std::pair<double, double> v1 = { p1,p3 };
       // std::pair<double, double> v2 = { p2,p4 };
       // std::pair<double, double> v3 = { p5,p7 };
       // std::pair<double, double> v4 = { p6,p8 };

       // points.push_back(krs::Vector2<double>{ p1, p3});
       // points.push_back(krs::Vector2<double>{ p2, p4});
       // points.push_back(krs::Vector2<double>{ p5, p7});
       // points.push_back(krs::Vector2<double>{ p6, p8});

       // medges.push_back({ v1, v3 });
       // medges.push_back({ v2, v4 });

       // std::cout << "contains working? " << krs::containsonePoint(medges.back(), v2, v1) << std::endl;
       // std::cout << "intersect working: " << delaunay.intersect(v1,v2,v3,v4) << std::endl;


    


    sf::RenderWindow window(sf::VideoMode(800, 600), "Delaunay triangulation");
    window.setFramerateLimit(1);

    // Transform each points of each vector as a rectangle
    for (const auto p : points) {
        sf::RectangleShape s{ sf::Vector2f(4, 4) };
        s.setPosition(static_cast<double>(p.x), static_cast<double>(p.y));
        window.draw(s);
    }

    

   for (const auto& e : edges) {
        sf::Vertex vertices[2] =
        {
            sf::Vertex(sf::Vector2f(e.first.first, e.first.second)),
            sf::Vertex(sf::Vector2f(e.second.first, e.second.second))
        };
        window.draw(vertices, 2, sf::Lines);
    }
   
    window.display();

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
    }
    return 0;
}