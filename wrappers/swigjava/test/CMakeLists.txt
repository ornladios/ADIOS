find_package(Java REQUIRED)
include(UseJava)

set (CMAKE_JAVA_TARGET_VERSION 1.8)

set (JAVA_TEST_SRC
  helloadios.java)

message (STATUS adiosjava_dir $<TARGET_FILE:adiosjava>)

get_target_property(ornl-adiosjava_jarfile ornl-adiosjava JAR_FILE)
get_target_property(ornl-adiosjava_classdir ornl-adiosjava CLASSDIR)

set(CMAKE_JAVA_INCLUDE_PATH ${ornl-adiosjava_jarfile})
add_jar (test-adiosjava
  ${JAVA_TEST_SRC})

get_target_property(test-adiosjava_jarfile test-adiosjava JAR_FILE)
get_target_property(test-adiosjava_classdir test-adiosjava CLASSDIR)

message(STATUS "Jar file ${adiosjava_jarfile}")
message(STATUS "Class compiled to ${adiosjava_classdir}")

enable_testing()
add_test(NAME test1 COMMAND ${Java_JAVA_EXECUTABLE} -Djava.library.path=$<TARGET_FILE_DIR:adiosjava> -cp ${ornl-adiosjava_jarfile}:${test-adiosjava_jarfile} helloadios)
