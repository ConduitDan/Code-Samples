#include <memory>
#include <cmath>
template<class Item>
class GridNeighborIterator;
template<class Item>
class SquareGridNeighborIterator;
template<class Item>
class TriangularGridNeighborIterator;


/* 
    A grid spatially organizes the items it holds into a lattice. Each element
    has an associated point in space and is connected to other points which can
    be accessed though a neighbor iterator. Different grid classes each can make
    its own nearest neighbor iterator and print out coordinates. 
*/
template<class Item>
class Grid {
protected:
    int _Nx;
    int _Ny;
    std::unique_ptr<Item[]> _values;
public:
    Grid(int Nx, int Ny):_Nx{Nx},_Ny{Ny}{
        _values = std::unique_ptr<Item[]>(new Item[get_num_elements()]);
    }

    virtual GridNeighborIterator<Item>* get_neighbor_iterator() = 0;
    Item& operator[](int index){return _values[index];}
    int get_num_elements(){return _Nx*_Ny;}
    virtual std::array<double,2> get_coord(int index) = 0;

    int get_Nx(){return _Nx;}
    int get_Ny(){return _Ny;}
};

template<class Item>
class SquareGrid: public Grid<Item> {
public:
    using Grid<Item>::Grid; // use only the constructor(s) from Grid
    GridNeighborIterator<Item>* get_neighbor_iterator() {
        return new SquareGridNeighborIterator<Item>(this);
    }
    std::array<double,2> get_coord(int index){
        int row_num = index/ this->_Nx;
        int col_num = index% this->_Nx;
        return {(double)col_num,(double)row_num};
    }
};

template<class Item>
class TriangularGrid: public Grid<Item> {
public:
    using Grid<Item>::Grid; // use only the constructor(s) from Grid
    GridNeighborIterator<Item>* get_neighbor_iterator(){
        return new TriangularGridNeighborIterator<Item>(this);
    }
    std::array<double,2> get_coord(int index){
        int row_num = index / this->_Nx;
        int col_num = index % this->_Nx;
        bool offset = (row_num)%2;
        return {col_num+0.5*offset,row_num*sqrt(3)/2};
    }
};







template<class Item>
class GridNeighborIterator {
protected:
    Grid<Item> * _parentGrid;
    int _current;
    int _start;
    int _numNeighbors;
    int _count;



public:
    GridNeighborIterator(Grid<Item> *parentGrid):_parentGrid{parentGrid}{}
    Item& operator*(void){return (*_parentGrid)[_current];}

    int get_index(){return _current;}

    virtual void calc_neighbor() = 0;

    virtual void begin(int index){_start = index; _count = 0; calc_neighbor();}
    void next() {_count++; calc_neighbor();}
    bool is_done() {return _count==_numNeighbors;}

};
template<class Item>
class SquareGridNeighborIterator: public GridNeighborIterator<Item> {
protected:
    using::GridNeighborIterator<Item>::_current;
    using::GridNeighborIterator<Item>::_numNeighbors;
    using::GridNeighborIterator<Item>::_start;
    using::GridNeighborIterator<Item>::_count;
    using::GridNeighborIterator<Item>::_parentGrid;
    
    int _neighborList[4];
    

public:
    SquareGridNeighborIterator(Grid<Item> *parentGrid):GridNeighborIterator<Item>(parentGrid){
        _numNeighbors = 4; 
    }
    void PeriodicBCLeft(int index){
        if (index % _parentGrid->get_Nx() == 0) {
            // if we're on the left wall add Nx to left
            _neighborList[0]+=_parentGrid->get_Nx();
        }
    }
    void PeriodicBCRight(int index){
        if ((index+1) % _parentGrid->get_Nx() == 0) {
            // if we're on the right wall subtract Nx to everything looking right
            _neighborList[1]-=_parentGrid->get_Nx();
        }
    }
    void PeriodicBCUp(int index){
        // if we're on the top subtract the number of total elements to everything looking up
        if (index/_parentGrid->get_Nx() == _parentGrid->get_Ny()-1){
            _neighborList[3]-=_parentGrid->get_num_elements();
        }
    }
    void PeriodicBCDown(int index){
        // if we're on the bottom add the number of total elements to everything looking down
        if (index/_parentGrid->get_Nx() == 0){
            _neighborList[2]+=_parentGrid->get_num_elements();
        }
    }

    void begin(int index){
        _neighborList[0]  = -1; // look left
        _neighborList[1] = 1; //look right
        _neighborList[2]  = -_parentGrid->get_Nx(); // look down
        _neighborList[3] = _parentGrid->get_Nx(); //look up

        // apply boundary conditions
        PeriodicBCLeft(index);
        PeriodicBCRight(index);
        PeriodicBCUp(index);
        PeriodicBCDown(index);

        GridNeighborIterator<Item>::begin(index); // call the parent function

    }
    void calc_neighbor(){
        if (_count<_numNeighbors){
            _current = _neighborList[_count]+_start;
        }
    }
};
template<class Item>
class TriangularGridNeighborIterator: public GridNeighborIterator<Item> {

protected:
    using::GridNeighborIterator<Item>::_current;
    using::GridNeighborIterator<Item>::_numNeighbors;
    using::GridNeighborIterator<Item>::_start;
    using::GridNeighborIterator<Item>::_count;
    using::GridNeighborIterator<Item>::_parentGrid;


    int _neighborList[6];
    bool _offsetRow; 
/*  +0  +1   2   3   3Nx-1
2Nx  *---*---*---*---  
      \0/ \1/ \2/ \ 2Nx-1
Nx +0  *---*---*---*-  offSet row looks up/down and left by +/- Nx
      / \ / \ / \ /                                right by +/- Nx +1
     *---*---*---*--- non offSet row look up/down and left by +/-Nx -1
    0    1    2   3                                  right by +/-Nx*/
public:

    TriangularGridNeighborIterator(Grid<Item> *parentGrid):GridNeighborIterator<Item>(parentGrid){
        _numNeighbors = 6;
    }
    void PeriodicBCLeft(int index){
        if (index % _parentGrid->get_Nx() == 0) {
            // if we're on the left wall add Nx to left
            _neighborList[0]+=_parentGrid->get_Nx();
            
            // if we're not off set we have to add to the up and down look left too
            if (!_offsetRow){
                _neighborList[2]+=_parentGrid->get_Nx();
                _neighborList[4]+=_parentGrid->get_Nx();
            }
        }
    }
    void PeriodicBCRight(int index){
        if ((index+1) % _parentGrid->get_Nx() == 0) {
            // if we're on the right wall subtract Nx to everything looking right
            _neighborList[1]-=_parentGrid->get_Nx();
            
            // if we're off set we have to subtract from the up and down look right too
            if (_offsetRow){
                _neighborList[3]-=_parentGrid->get_Nx();
                _neighborList[5]-=_parentGrid->get_Nx();
            }
        }
    }
    void PeriodicBCUp(int index){
        // if we're on the top subtract the number of total elements to everything looking up
        if (index/_parentGrid->get_Nx() == _parentGrid->get_Nx()-1){
            _neighborList[2]-=_parentGrid->get_num_elements();
            _neighborList[3]-=_parentGrid->get_num_elements();
        }
    }
    void PeriodicBCDown(int index){
        // if we're on the bottom add the number of total elements to everything looking down
        if (index/_parentGrid->get_Nx() == 0){
            _neighborList[4]+=_parentGrid->get_num_elements();
            _neighborList[5]+=_parentGrid->get_num_elements();
        }
    }
    void begin(int index){
        // every other row is an offset row
        // _start/_parentGrid->get_Nx() gets the row number
        _neighborList[0]  = -1; // look left
        _neighborList[1] = 1; //look right
        _neighborList[2]  = _parentGrid->get_Nx() - 1; // look up left (non offset)
        _neighborList[3] = _parentGrid->get_Nx(); //look up right (non offset)
        _neighborList[4]  = -_parentGrid->get_Nx() -1; // look down left (non offset)
        _neighborList[5] = -_parentGrid->get_Nx(); //look down right (non offset)

        _offsetRow = (index/_parentGrid->get_Nx())%2;
        // if we're on an off set row ad one to the look up and down entires
        if (_offsetRow){
            for (int i = 2; i<6; i++){
                _neighborList[i]++;
            }

        }
        PeriodicBCLeft(index);
        PeriodicBCRight(index);
        PeriodicBCUp(index);
        PeriodicBCDown(index);

        GridNeighborIterator<Item>::begin(index); // call the parent function

    }
    void calc_neighbor(){
        if (_count<_numNeighbors){
            _current =_start + _neighborList[_count];
         }
    }
};

