// $Id: t_db.cc,v 1.2 2008-07-30 19:27:11 edwards Exp $
/*! \file
 *  \brief Test the database routines
 */

#include "chroma.h"
#include "ConfVarDSizeStoreDB.h"

namespace Chroma
{
  using namespace FFDB;

  // Simple concrete key class
  class TestDBKey : public DBKey
  {
  public:
    // Part of Serializable
    const unsigned short serialID (void) const {return 789;}
    void objectBuffer (char* &buffer, unsigned int& length) {
      buffer = (char*)data_; length = sizeof(data_);
    }

    const unsigned int sizeOfObject (void) const throw (SerializeException) {
      return sizeof(data_);
    }

    void writeObject (void) throw (SerializeException) {
      data_[0] = 1;
      data_[1] = 17;
    }

    void readObject (void) throw (SerializeException) {
    }

    // Part of DBKey
    int hasHashFunc (void) const {return 0;}
    int hasCompareFunc (void) const {return 0;}

    /**
     * Static Hash Function Implementation
     */
    static unsigned int hash (Db *db, const void* bytes, unsigned int len) {return 0;}

    /**
     * Static empty compare function 
     */
    static int compare (Db *db, const Dbt* k1, const Dbt* k2) {return 0;}

  private:
    int a;
    int b;
    int data_[2];
  };


  // Simple concrete data class
  class TestDBData : public DBData
  {
  public:
    const unsigned short serialID (void) const {return 890;}
    void objectBuffer (char* &buffer, unsigned int& length) {
      buffer = (char*)data_; length = sizeOfObject();
    }

    const unsigned int sizeOfObject (void) const throw (SerializeException) {
      return sizeof(data_);
    }

    void writeObject (void) throw (SerializeException) {
      for(int i=0; i < 10; ++i)
	data_[i] =  i;
    }

    void readObject (void) throw (SerializeException) {
    }


  private:
    float data_[10];
  };

}


using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_db.xml");
  push(xml,"t_db");

//  push(xml,"lattice");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
//  pop(xml);

  // Try out some simple DB stuff
  try
  {
    TestDBKey testDBKey;
    TestDBData testDBData;

    // Open it
    QDPIO::cout << "open" << endl;
    ConfVarDSizeStoreDB<TestDBKey, TestDBData> db("test.db");

    // Test it
    QDPIO::cout << "insert" << endl;
    db.insert(testDBKey, testDBData);

    // Flush it
    QDPIO::cout << "flush" << endl;
    db.flush();

    QDPIO::cout << "closing" << endl;
  }
  catch (...) {
    QDPIO::cout << "Error in db routines" << endl;
  }

  // Test the db
  try
  {
    // Open it
    QDPIO::cout << "open" << endl;
    typedef ConfVarDSizeStoreDB<TestDBKey, TestDBData> DBType_t;
    DBType_t db("test.db");

    // Test it
    QDPIO::cout << "find keys" << endl;
    std::vector<TestDBKey> keys;
    db.keys(keys);

    // Flush it
    QDPIO::cout << "printing" << endl;
    for(int i=0; i < keys.size(); ++i)
      QDPIO::cout << "found a key: i= " << i << endl;
      
//    for(DBType_t::key_iterator ptr=db.begin(); ptr != db.end(); ++ptr)
//    {
//      QDPIO::cout << "found a key" << endl;
//    }

    QDPIO::cout << "closing" << endl;
  }
  catch (DbException& e) {
    QDPIO::cerr << "DBException: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch(std::exception &e) {
    QDPIO::cerr << "Std exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch (...) {
    QDPIO::cout << "Generic error in db routines" << endl;
    QDP_abort(1);
  }

  // Done
  pop(xml);

  // Time to bolt
  Chroma::finalize();

  QDPIO::cout << "finished" << endl;
  exit(0);
}

