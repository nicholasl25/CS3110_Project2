open OUnit2
open Helices
open Dna

(******************************************************************************)
(********************** REMEMBER TO WRITE TESTS FIRST! ************************)
(******************************************************************************)

(* Pretty printers for OUnit messages. These print values mostly like the
   toplevel prints them. *)

let pp_pair pp_l pp_r (l, r) = Printf.sprintf "(%s, %s)" (pp_l l) (pp_r r)

let pp_list pp_elt lst =
  let pp_elts lst =
    let rec loop n acc = function
      | [] -> acc
      | [ h ] -> acc ^ pp_elt h
      | h1 :: (h2 :: t as t') ->
          if n = 100 then acc ^ "..." (* stop printing long list *)
          else loop (n + 1) (acc ^ pp_elt h1 ^ "; ") t'
    in
    loop 0 "" lst
  in
  Printf.sprintf "[%s]" (pp_elts lst)

let pp_helix = pp_list string_of_nucleotide
let pp_helix_list = pp_list pp_helix
let pp_acid_list = pp_list string_of_acid
let pp_acid_list_list = pp_list pp_acid_list

let pp_tree t =
  let rec f ws t =
    let ews = ws ^ "    " in
    match t with
    | Leaf h -> ws ^ "Leaf " ^ pp_helix h
    | Node (l, r) ->
        ws ^ "Node (\n" ^ f ews l ^ ",\n" ^ f ews r ^ "\n" ^ ws ^ ")"
  in
  "\n" ^ f "" t

let pp_ltree t =
  let rec f ws t =
    let ews = ws ^ "    " in
    match t with
    | LLeaf h -> ws ^ "LLeaf " ^ pp_helix h
    | LNode (l, h, r) ->
        ws ^ "LNode (\n" ^ f ews l ^ ",\n" ^ ews ^ pp_helix h ^ ",\n" ^ f ews r
        ^ "\n" ^ ws ^ ")"
  in
  "\n" ^ f "" t

let pp_ltree_list = pp_list pp_ltree
let pp_ltree_int_pair = pp_pair pp_ltree string_of_int

(************************ complementary_helix tests ***************************)

let complementary_helix_test out in1 _ =
  assert_equal ~printer:pp_helix
    ~msg:("function: complementary_helix\ninput: " ^ pp_helix in1)
    out
    (Dna.complementary_helix in1)

let complementary_helix_tests =
  [
    "complementary_helix multi-element list"
    >:: complementary_helix_test [ G; A; A; G ] [ C; T; T; C ];
  ]

(*************************** hamming_distance tests ***************************)

let hamming_distance_test out in1 in2 _ =
  assert_equal ~printer:string_of_int
    ~msg:
      (Printf.sprintf "function: hamming_distance\ninputs: %s %s" (pp_helix in1)
         (pp_helix in2))
    out
    (Dna.hamming_distance in1 in2)

let hamming_invalid_arg_test in1 in2 _ =
  let exn = Invalid_argument "..." in
  assert_raises
    ~msg:
      (Printf.sprintf "function: hamming_distance\ninput: %s %s" (pp_helix in1)
         (pp_helix in2))
    exn
    (fun () ->
      try Dna.hamming_distance in1 in2 with Invalid_argument _ -> raise exn)

let hamming_distance_tests =
  [ (* Example test cases: *)
    (* "hamming_distance empty" >:: hamming_distance_test 0 [] []; *)
    (* "hamming_distance one different nucleotide" >:: hamming_distance_test 1 [
       G; A ] [ G; T ]; *)
    (* "hamming_distance nonmatching lengths" >:: hamming_invalid_arg_test [ C ]
       []; *) ]

(************************ decreasing_similarity tests *************************)

let decreasing_similarity_test out in1 _ =
  assert_equal ~printer:string_of_bool
    ~msg:("function: decreasing_similarity\ninput: " ^ pp_helix_list in1)
    out
    (Dna.decreasing_similarity in1)

let decreasing_similarity_tests =
  [ (* Example test case: *)
    (* "decreasing_similarity to humans" >:: decreasing_similarity_test true
       (Dna.most_like_human ()); *) ]

(***************************** count_leaves tests *****************************)

let count_leaves_test out in1 _ =
  assert_equal ~printer:string_of_int
    ~msg:("function: count_leaves\ninput: " ^ pp_tree in1)
    out (Dna.count_leaves in1)

let count_leaves_tests = []

(**************************** helix_of_tree tests *****************************)

let helix_of_tree_test out in1 _ =
  assert_equal ~printer:pp_helix
    ~msg:("function: helix_of_tree\ninput: " ^ pp_ltree in1)
    out (Dna.helix_of_tree in1)

let helix_of_tree_tests = []

(**************************** unlabel_tree tests ******************************)

let unlabel_tree_test out in1 _ =
  assert_equal ~printer:pp_tree
    ~msg:("function: unlabel_tree\ninput: " ^ pp_ltree in1)
    out (Dna.unlabel_tree in1)

let unlabel_tree_tests = []

(************************* guess_parent_helix tests ***************************)

let guess_parent_helix_test out in1 in2 _ =
  assert_equal ~printer:pp_helix
    ~msg:
      (Printf.sprintf "function: guess_parent_helix\ninputs: %s %s"
         (pp_helix in1) (pp_helix in2))
    out
    (Dna.guess_parent_helix in1 in2)

let guess_invalid_arg_test in1 in2 _ =
  let exn = Invalid_argument "..." in
  assert_raises
    ~msg:
      (Printf.sprintf "function: guess_parent_helix\ninput: %s %s"
         (pp_helix in1) (pp_helix in2))
    exn
    (fun () ->
      try Dna.guess_parent_helix in1 in2 with Invalid_argument _ -> raise exn)

let guess_parent_helix_tests =
  [ (* Example test case: *)
    (* "guess_parent_helix one difference" >:: guess_parent_helix_test [ G; C; A
       ] [ T; C; A ] [ G; C; A ]; *) ]

(************************ add_ancestor_labels tests ***************************)

let add_ancestor_labels_test out in1 _ =
  assert_equal ~printer:pp_ltree
    ~msg:("function: add_ancestor_labels\ninput: " ^ pp_tree in1)
    out
    (Dna.add_ancestor_labels in1)

let add_ancestor_labels_tests =
  [ (* Example test case: *)
    (* "add_ancestor_labels leaf" >:: add_ancestor_labels_test (LLeaf [ T; C ])
       (Leaf [ T; C ]); *) ]

(************************ parent_child_hamming tests **************************)

(* Be sure to test for trees of depth greater than one. *)

let parent_child_hamming_test out in1 _ =
  assert_equal ~printer:string_of_int
    ~msg:("function: parent_child_hamming\ninput: " ^ pp_ltree in1)
    out
    (Dna.parent_child_hamming in1)

let parent_child_hamming_tests =
  [ (* Example test case: *)
    (* "parent_child_hamming depth-2 tree, all different" >::
       parent_child_hamming_test 2 (LNode (LLeaf [ T ], [ A ], LLeaf [ G ])); *) ]

(**************************** simplest_tree tests *****************************)

let t1 = Dna.LNode (Dna.LLeaf [ A ], [ T ], Dna.LLeaf [ C ])

let t2 =
  Dna.LNode
    ( Dna.LLeaf [ G; T ],
      [ A; T ],
      Dna.LNode (Dna.LLeaf [ T; T ], [ T; T ], Dna.LLeaf [ C; G ]) )

let simplest_tree_test out in1 _ =
  assert_equal ~printer:pp_ltree_int_pair
    ~msg:("function: simplest_tree\ninput: " ^ pp_ltree_list in1)
    out (Dna.simplest_tree in1)

let simplest_tree_invalid_arg_test in1 _ =
  let exn = Invalid_argument "..." in
  assert_raises
    ~msg:(Printf.sprintf "function: simplest_tree\ninput: " ^ pp_ltree_list in1)
    exn
    (fun () -> try Dna.simplest_tree in1 with Invalid_argument _ -> raise exn)

let simplest_tree_tests =
  [ (* Example test case: *)
    (* "simplest_tree two tree list" >:: simplest_tree_test (t1, 2) [ t1; t2 ]  *) ]

(********************* find_simplest_tree tests ************************)

let find_simplest_tree_test out in1 in2 in3 in4 _ =
  assert_equal ~printer:pp_tree
    ~msg:
      (Printf.sprintf "function: find_simplest_tree\ninputs: %s %s %s %s"
         (pp_helix in1) (pp_helix in2) (pp_helix in3) (pp_helix in4))
    out
    (Dna.find_simplest_tree in1 in2 in3 in4)

let find_simplest_invalid_arg_test in1 in2 in3 in4 _ =
  let exn = Invalid_argument "..." in
  assert_raises
    ~msg:
      (Printf.sprintf "function: find_simplest_tree\ninputs: %s %s %s %s"
         (pp_helix in1) (pp_helix in2) (pp_helix in3) (pp_helix in4))
    exn
    (fun () ->
      try Dna.find_simplest_tree in1 in2 in3 in4
      with Invalid_argument _ -> raise exn)

let find_simplest_tree_tests =
  [ (* Example test case: *)
    (* "simplest_greater_ape_tree" >:: find_simplest_tree_test
       (Dna.simplest_greater_ape_tree ()) gorilla human chimpanzee orangutan; *) ]

(**************************** acids_of_helix tests ****************************)

let acids_of_helix_test out in1 _ =
  assert_equal ~printer:pp_acid_list
    ~msg:("function: acids_of_helix\ninput: " ^ pp_helix in1)
    out (Dna.acids_of_helix in1)

let acids_of_helix_tests =
  [ (* Example test case: *)
    (* "acids_of_helix single codon" >:: acids_of_helix_test [ Met ] [ A; G; T;
       A; C ]; *) ]

(*************************** all_acids_of_helix tests *************************)

let all_acids_of_helix_test out in1 _ =
  assert_equal ~printer:pp_acid_list_list
    ~msg:("function: all_acids_of_helix\ninput: " ^ pp_helix in1)
    out
    (Dna.all_acids_of_helix in1)

let all_acids_of_helix_tests =
  [ (* Example test case: *)
    (* "all_acids_of_helix [T; A; C; A; C; T]" >:: all_acids_of_helix_test [ [
       Met ] ] [ T; A; C; A; C; T ]; *) ]

(******************************************************************************)

let tests =
  "dna test suite"
  >::: complementary_helix_tests @ hamming_distance_tests
       @ decreasing_similarity_tests @ count_leaves_tests @ helix_of_tree_tests
       @ unlabel_tree_tests @ guess_parent_helix_tests
       @ add_ancestor_labels_tests @ parent_child_hamming_tests
       @ simplest_tree_tests @ find_simplest_tree_tests @ acids_of_helix_tests
       @ all_acids_of_helix_tests

let _ = run_test_tt_main tests
