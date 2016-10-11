/*
 * char **STRING_ARRAY typemaps.
 * These typemaps are for C String arrays which are NULL terminated.
 *   char *values[] = { "one", "two", "three", NULL }; // note NULL
 * char ** is mapped to a Java String[].
 *
 * Example usage wrapping:
 *   %apply char **STRING_ARRAY { char **input };
 *   char ** foo(char **input);
 *
 * Java usage:
 *   String numbers[] = { "one", "two", "three" };
 *   String[] ret = modulename.foo( numbers };
 */
%typemap(jni) char **STRING_ARRAY "jobjectArray"
%typemap(jtype) char **STRING_ARRAY "String[]"
%typemap(jstype) char **STRING_ARRAY "String[]"
%typemap(in) char **STRING_ARRAY (jint size) {
  int i = 0;
  if ($input) {
    size = JCALL1(GetArrayLength, jenv, $input);
#ifdef __cplusplus
    $1 = new char*[size+1];
#else
    $1 = (char **)malloc((size+1) * sizeof(char *));
#endif
    for (i = 0; i<size; i++) {
      jstring j_string = (jstring)JCALL2(GetObjectArrayElement, jenv, $input, i);
      const char *c_string = JCALL2(GetStringUTFChars, jenv, j_string, 0);
#ifdef __cplusplus
      $1[i] = new char [strlen(c_string)+1];
#else
      $1[i] = (char *)malloc((strlen(c_string)+1) * sizeof(const char *));
#endif
      strcpy($1[i], c_string);
      JCALL2(ReleaseStringUTFChars, jenv, j_string, c_string);
      JCALL1(DeleteLocalRef, jenv, j_string);
    }
    $1[i] = 0;
  } else {
    $1 = 0;
    size = 0;
  }
}

%typemap(freearg) char **STRING_ARRAY {
  int i;
  for (i=0; i<size$argnum; i++)
#ifdef __cplusplus
    delete[] $1[i];
  delete[] $1;
#else
  free($1[i]);
  free($1);
#endif
}

%typemap(out) char **STRING_ARRAY {
  if ($1) {
    int i;
    jsize len=0;
    jstring temp_string;
    const jclass clazz = JCALL1(FindClass, jenv, "java/lang/String");

    //while ($1[len]) len++;
    len=arg1->nvars;
    $result = JCALL3(NewObjectArray, jenv, len, clazz, NULL);
    /* exception checking omitted */

    for (i=0; i<len; i++) {
      //printf("%d, %s\n", i, $1[i]);
      temp_string = JCALL1(NewStringUTF, jenv, $1[i]);
      JCALL3(SetObjectArrayElement, jenv, $result, i, temp_string);
      JCALL1(DeleteLocalRef, jenv, temp_string);
    }
  }
}

%typemap(javain) char **STRING_ARRAY "$javainput"
%typemap(javaout) char **STRING_ARRAY {
    return $jnicall;
}
