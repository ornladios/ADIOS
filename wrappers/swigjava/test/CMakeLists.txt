find_package(Java REQUIRED)
include(UseJava)

set (CMAKE_JAVA_TARGET_VERSION 1.8)

set (JAVA_TEST_SRC
  helloadioswriter.java
  helloadios.java)

message (STATUS adiosjava_dir $<TARGET_FILE:adiosjava>)

get_target_property(ornl-adiosjava_jarfile ornl-adiosjava JAR_FILE)
get_target_property(ornl-adiosjava_classdir ornl-adiosjava CLASSDIR)

set(CMAKE_JAVA_INCLUDE_PATH ${ornl-adiosjava_jarfile})
add_jar (test-adiosjava
  SOURCES ${JAVA_TEST_SRC})

get_target_property(test-adiosjava_jarfile test-adiosjava JAR_FILE)
get_target_property(test-adiosjava_classdir test-adiosjava CLASSDIR)

message(STATUS "Jar file ${adiosjava_jarfile}")
message(STATUS "Class compiled to ${adiosjava_classdir}")

configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/arrays.xml
  ${CMAKE_CURRENT_BINARY_DIR}/arrays.xml
  )

enable_testing()
add_test(NAME test1 COMMAND ${Java_JAVA_EXECUTABLE} -Djava.library.path=$<TARGET_FILE_DIR:adiosjava> -cp ${ornl-adiosjava_jarfile}:${test-adiosjava_jarfile} helloadioswriter)
add_test(NAME test2 COMMAND ${Java_JAVA_EXECUTABLE} -Djava.library.path=$<TARGET_FILE_DIR:adiosjava> -cp ${ornl-adiosjava_jarfile}:${test-adiosjava_jarfile} helloadios)
