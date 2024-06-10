#include "Grid.h"
#define BOOST_TEST_MODULE "Grid Tests"
#include <boost/test/included/unit_test.hpp>
#include <set>

struct NeighborFixture {
    std::set<int> correct;
    std::set<int> code_answer;
    TriangularGrid<double> g;
    NeighborFixture():g(4,4){}
    void test(int index) {
        auto it = g.get_neighbor_iterator();
        for (it->begin(index); !it->is_done(); it->next()){
            code_answer.insert(it->get_index());
        }
        BOOST_CHECK(code_answer==correct);
        if (code_answer!=correct){
            for (auto node : code_answer) std::cout<<node<<' ';
            std::cout<<std::endl;
            std::cout<<"Expected: "<<std::endl;
            for (auto node : correct) std::cout<<node<<' ';
            std::cout<<std::endl;

        }
        if (code_answer!=correct){
        }

    }
};

BOOST_AUTO_TEST_SUITE(TriangularLatticeIterator);


BOOST_FIXTURE_TEST_CASE(middle_offset,NeighborFixture){
    correct = {4,6,9,10,1,2};
    test(5);

}
BOOST_FIXTURE_TEST_CASE(middle,NeighborFixture){
    correct = {8,4,5,10,13,12};
    test(9);
}

BOOST_FIXTURE_TEST_CASE(left_bottom_wall,NeighborFixture){
    correct = {1, 4, 3, 7, 15, 12};
    test(0);
}
BOOST_FIXTURE_TEST_CASE(right_bottom_wall,NeighborFixture){
    correct = {2, 6, 7, 0, 15, 14};
    test(3);
}
BOOST_FIXTURE_TEST_CASE(left_offset_wall,NeighborFixture){
    correct = {0, 1, 5, 8, 9, 7};
    test(4);
}
BOOST_FIXTURE_TEST_CASE(right_offet_wall,NeighborFixture){
    correct = {3, 6, 11, 0, 8, 4};
    test(7);
}
BOOST_FIXTURE_TEST_CASE(left_top_wall,NeighborFixture){
    correct = {8, 9, 13, 0, 15, 1};
    test(12);
}
BOOST_FIXTURE_TEST_CASE(right_top_wall,NeighborFixture){
    correct = {14, 12, 11, 8, 3, 0};
    test(15);
}

BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE(SquareLattice);
BOOST_AUTO_TEST_CASE(right_top_wall){
    SquareGrid<double> g(4,4);
    auto it = g.get_neighbor_iterator();

    
    std::set<int> correct = {3,12,14,11};
    std::set<int> code_answer;
    for (it->begin(15); !it->is_done(); it->next()){
        code_answer.insert(it->get_index());
    }
    BOOST_CHECK(code_answer==correct);
    if (code_answer!=correct){
        for (auto node : code_answer) std::cout<<node<<' ';
        std::cout<<std::endl;
    }
}
BOOST_AUTO_TEST_CASE(right_bottom_wall){
    SquareGrid<double> g(4,4);
    auto it = g.get_neighbor_iterator();

    
    std::set<int> correct = {7,2,0,15};
    std::set<int> code_answer;
    for (it->begin(3); !it->is_done(); it->next()){
        code_answer.insert(it->get_index());
    }
    BOOST_CHECK(code_answer==correct);
    if (code_answer!=correct){
        for (auto node : code_answer) std::cout<<node<<' ';
        std::cout<<std::endl;
    }


}
BOOST_AUTO_TEST_CASE(left_top_wall){
    SquareGrid<double> g(4,4);
    auto it = g.get_neighbor_iterator();

    
    std::set<int> correct = {13,8,0,15};
    std::set<int> code_answer;
    for (it->begin(12); !it->is_done(); it->next()){
        code_answer.insert(it->get_index());
    }
    BOOST_CHECK(code_answer==correct);
    if (code_answer!=correct){
        for (auto node : code_answer) std::cout<<node<<' ';
        std::cout<<std::endl;
    }


}
BOOST_AUTO_TEST_CASE(left_bottom_wall){
    SquareGrid<double> g(4,4);
    auto it = g.get_neighbor_iterator();

    
    std::set<int> correct = {1,4,3,12};
    std::set<int> code_answer;
    for (it->begin(0); !it->is_done(); it->next()){
        code_answer.insert(it->get_index());
    }
    BOOST_CHECK(code_answer==correct);
    if (code_answer!=correct){
        for (auto node : code_answer) std::cout<<node<<' ';
        std::cout<<std::endl;
    }


}

BOOST_AUTO_TEST_SUITE_END();
