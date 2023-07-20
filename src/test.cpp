#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE VALIDATE_P452
#include "boost/test/unit_test.hpp"
#include "P452.hpp"
#include "P452.test.hpp"

BOOST_AUTO_TEST_CASE( p452_test1 )
{
    BOOST_CHECK( test1() );        // #1 continues on error
}

BOOST_AUTO_TEST_CASE( p452_test2 )
{
    BOOST_CHECK( test2() );        // #1 continues on error
}

BOOST_AUTO_TEST_CASE( p452_test3 )
{
    BOOST_CHECK( test3() );        // #1 continues on error
}

BOOST_AUTO_TEST_CASE( p452_test4 )
{
    BOOST_CHECK( test4() );        // #1 continues on error
}

BOOST_AUTO_TEST_CASE( p452_test5 )
{
    BOOST_CHECK( test5() );        // #1 continues on error
}

BOOST_AUTO_TEST_CASE( p452_test6 )
{
    BOOST_CHECK( test6() );        // #1 continues on error
}

BOOST_AUTO_TEST_CASE( p452_test7 )
{
    BOOST_CHECK( test7() );        // #1 continues on error
}

BOOST_AUTO_TEST_CASE( my_test8 )
{

    BOOST_CHECK( test8() );        // #1 continues on error

}
// testing reading and interpolation from the meteorological maps
BOOST_AUTO_TEST_CASE( my_test9 )
{

    BOOST_CHECK( test9() );        // #1 continues on error

}

// testing great circle path computations (mid path)
BOOST_AUTO_TEST_CASE( my_test10 )
{
    BOOST_CHECK( test10() );        // #1 continues on error

}

// testing great circle path computations (mid path)
BOOST_AUTO_TEST_CASE( my_test11 )
{
    BOOST_CHECK( test11() );        // #1 continues on error

}

// testing great circle path computations (mid path)
BOOST_AUTO_TEST_CASE( my_test12 )
{
    BOOST_CHECK( test12() );        // #1 continues on error

}
