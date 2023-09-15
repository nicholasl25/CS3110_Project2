open Helices

(** A list of ape species and their corresponding names. *)
let ape_names : (helix * string) list =
  [
    (lar_gibbon, "Lar Gibbon");
    (pileated_gibbon, "Pileated Gibbon");
    (siamang, "Siamang");
    (white_cheeked_gibbon, "White Cheeked Gibbon");
    (orangutan, "Orangutan");
    (gorilla, "Gorilla");
    (chimpanzee, "Chimpanzee");
    (human, "Human");
  ]

(** Given a helix, returns the name of the ape it corresponds to (e.g.
    "Orangutan") as given in ape_names. Returns "Not an Ape!" if the helix does
    not match any of the ape helices in ape_names. *)
let string_of_ape (ape : helix) : string =
  let rec lookup_ape (l : (helix * string) list) : string =
    match l with
    | (h, name) :: t -> if ape = h then name else lookup_ape t
    | [] -> "Not an Ape!"
  in
  lookup_ape ape_names

(** Returns a string representation of a given nucleotide. *)
let string_of_nucleotide = function
  | A -> "A"
  | T -> "T"
  | G -> "G"
  | C -> "C"

(** Given a single helix, computes the complementary helix of the double helix.
    The complementary helix is obtained by swapping A <-> T and G <-> C. *)
let swap (x: nucleotide) : nucleotide =
  match x with
  | A -> T
  | T -> A
  | G -> C
  | C -> G
    
let rec complementary_helix (x : helix) : helix = 
  match x with
    |[] -> []
    |[h] -> [swap h]
    | h :: t -> swap h :: complementary_helix t
(** Given two equal-length helices, computes the Hamming distance between them.
    This is the number of corresponding positions in the helices at which the
    nucleotides in those positions differ. Raises: Invalid_argument if the input
    helices are not the same length. *)
let rec hamming_distance (x1 : helix) (x2 : helix) : int =
  if List.length x1 = List.length x2
    then match x1, x2 with
  | [], [] -> 0
  | [h], [k]-> if h = k then 0 else 1
  | h :: t1 , k :: t2 -> if h = k then hamming_distance t1 t2
    else (hamming_distance t1 t2) + 1
  else invalid_arg "Two Helixs are not of the same length"
 

(** A list of helices of ape species in order of non-increasing similarity
    (highest to lowest) to humans, as measured by Hamming distance. You can
    break ties however you want. *)
let most_like_human () : helix list = 
  [chimpanzee; gorilla; orangutan; siamang; pileated_gibbon;white_cheeked_gibbon; lar_gibbon]

(** Determines whether a list of ape helices is ordered in decreasing order of
    similarity to humans, as measured by Hamming distance. Returns true for
    empty lists and lists that contain exactly one helix. *)

let human_distance (ape: helix) : int = 
  hamming_distance human ape
let rec decreasing_similarity (apes : helix list) : bool = 
  match apes with 
  |[] -> true
  |[h] -> true
  | h :: [t] -> human_distance h <= human_distance t
  | h :: k :: t -> human_distance h <= human_distance k && decreasing_similarity (k :: t)
  

(** An evolutionary tree is a binary tree whose leaves are labeled with DNA
    helices. *)
type tree =
  | Leaf of helix
  | Node of tree * tree

(* A representation of the left side of the tree from Figure A in the A1 writeup
   (the "greater apes"). *)
let greater_apes () : tree =
  Node (Leaf orangutan, Node (Leaf gorilla, Node (Leaf human, Leaf chimpanzee)))

(* A representation of the right side of the tree from Figure A in the A1
   writeup (the "lesser apes"). *)
let lesser_apes () : tree = 
  Node (Leaf white_cheeked_gibbon, Node (Leaf siamang, Node (Leaf lar_gibbon, Leaf pileated_gibbon)))

(** Counts the number of leaves of a tree. *)
let rec count_leaves (t : tree) : int = 
  match t with 
  |Leaf(_) -> 1
  |Node(l, r) -> count_leaves l + count_leaves r



(* Given four helices that should appear on the leaves of each output tree,
   returns a list of all distinct trees with these four input helices at their
   leaves. Two trees that are obtained from each other by permuting the children
   of any node are considered the same tree. *)
let rec all_trees (x1 : helix) (x2 : helix) (x3 : helix) (x4 : helix) :
    tree list =
  let tree1 (x1 : helix) (x2 : helix) (x3 : helix) (x4 : helix) : tree =
    Node (Node (Leaf x1, Leaf x2), Node (Leaf x3, Leaf x4))
  in
  let tree2 (x1 : helix) (x2 : helix) (x3 : helix) (x4 : helix) : tree =
    Node (Leaf x1, Node (Leaf x2, Node (Leaf x3, Leaf x4)))
  in
  [
    tree1 x1 x2 x3 x4;
    tree1 x1 x3 x2 x4;
    tree1 x1 x4 x2 x3;
    tree2 x1 x2 x3 x4;
    tree2 x1 x3 x2 x4;
    tree2 x1 x4 x2 x3;
    tree2 x2 x1 x3 x4;
    tree2 x2 x3 x1 x4;
    tree2 x2 x4 x1 x3;
    tree2 x3 x1 x2 x4;
    tree2 x3 x2 x1 x4;
    tree2 x3 x4 x1 x2;
    tree2 x4 x1 x2 x3;
    tree2 x4 x2 x1 x3;
    tree2 x4 x3 x1 x2;
  ]

(** Creates a list of all trees containing the four greater apes at the leaves. *)
let all_greater_ape_trees () : tree list = 
  all_trees orangutan gorilla human chimpanzee

(** Creates a list of all trees containing the four lesser apes at the leaves. *)
let all_lesser_ape_trees () : tree list = 
  all_trees white_cheeked_gibbon siamang lar_gibbon pileated_gibbon

(** A labeled_tree is a tree with a helix at each internal node representing the
    DNA of a possible ancestor. *)
type labeled_tree =
  | LLeaf of helix
  | LNode of labeled_tree * helix * labeled_tree

(** Returns the helix at the root node of a given labeled_tree. *)
let helix_of_tree (t : labeled_tree) : helix =
  match t with 
  |LLeaf(x) -> x
  |LNode(_,x,_) -> x

(** Given a labeled_tree, returns a tree of identical structure with the helix
    in each Node removed. *)
let rec unlabel_tree (t : labeled_tree) : tree = 
  match t with
  |LLeaf(x) -> Leaf(x)
  |LNode(l, _, r) -> Node(unlabel_tree l, unlabel_tree r)

  (** Helper functtion for guess_parent_helix, given two nucleotides this function 
    returns the first nucleotide if the nucleotides are equivalent and returns 
    the nucleotide that would occur first in dictionary order (A, C, G, T)
    if the nucleotides are different*)
  
let original_nucleotide (x1 : nucleotide) (x2 : nucleotide) =
  if x1 = x2 then x1 else
  match (x1, x2) with 
  | (A, _) -> A
  | (_, A) -> A
  | (C, _) -> C
  | (_, C) -> C
  | (_, _) -> G


(** Given two helices of equal length, returns a guess of a valid helix that
    could label their parent node in an evolutionary tree. The guessed helix is
    computed using nucleotide-by-nucleotide comparison of x1 and x2. If the
    nucleotides match, the corresponding parent nucleotide is the same. If they
    do not match, the corresponding parent nucleotide is assigned the one that
    would occur first in dictionary order (A, C, G, T). Raises: Invalid_argument
    if the given helices are not the same length. *)
let rec guess_parent_helix (x1 : helix) (x2 : helix) : helix =
  match (x1, x2) with 
  | ([], []) -> []
  | ([h1], [h2]) -> [original_nucleotide h1 h2]
  | (h1 :: t1, h2 :: t2) -> original_nucleotide h1 h2 :: guess_parent_helix t1 t2
  | (_, _) -> raise (Invalid_argument "Helixes are not of the same length")

(** Given an unlabeled tree, generates a labeled_tree with the same shape but
    with all internal nodes properly labeled with guessed helices at the
    internal nodes. Use guess_parent_helix to compute the appropriate label for
    each parent node, given the labels of its two children. *)
let rec add_ancestor_labels (t : tree) : labeled_tree = 
  match t with
  |Leaf(x) -> LLeaf(x)
  |Node(n1, n2) -> LNode(add_ancestor_labels n1, 
  guess_parent_helix (helix_of_tree (add_ancestor_labels n1)) (helix_of_tree (add_ancestor_labels n2)), 
  add_ancestor_labels n2)
     

(** Add labels to all trees in a list using add_ancestor_labels. *)
let rec add_ancestor_labels_list (ts : tree list) : labeled_tree list =
  match ts with
  | [] -> []
  | hd :: tl -> add_ancestor_labels hd :: add_ancestor_labels_list tl

(** Creates labeled_trees for lesser and greater apes. *)
let labeled_greater_ape_trees () : labeled_tree list =
  add_ancestor_labels_list (all_greater_ape_trees ())

let labeled_lesser_ape_trees () : labeled_tree list =
  add_ancestor_labels_list (all_lesser_ape_trees ())

(** Given a labeled_tree, computes the sum of Hamming distances between all
    parent-child pairs, where "parent-child pair" refers to a DIRECT
    parent-child relationship, not just an ancestor relationship. *)
let rec parent_child_hamming (t : labeled_tree) : int = 
  match t with
  |LLeaf(h) -> 0
  |LNode(l, h, r) -> hamming_distance h (helix_of_tree l) + hamming_distance h (helix_of_tree r) + 
  parent_child_hamming l + parent_child_hamming r

(** Given a list of labeled trees, returns the one in the list of lowest
    complexity along with its complexity, as measured by the
    parent_child_hamming function. If there are more than one tree with the same
    lowest complexity, return the first one in the list. Raises:
    Invalid_argument if the input list is empty. *)
let rec simplest_tree (ts : labeled_tree list) : labeled_tree * int =
  match ts with 
  |[] -> raise (Invalid_argument "Empty List")
  |[h] -> (h, parent_child_hamming h)
  |h::t -> if parent_child_hamming h <= snd(simplest_tree t) then (h, parent_child_hamming h)
  else simplest_tree t


(** Calculates the simplest evolutionary tree for the greater apes. *)
let simplest_greater_ape_tree () : tree =
  let t, _ = simplest_tree (labeled_greater_ape_trees ()) in
  unlabel_tree t

(** Calculates the simplest evolutionary tree for the lesser apes. *)
let simplest_lesser_ape_tree () : tree =
  let t, _ = simplest_tree (labeled_lesser_ape_trees ()) in
  unlabel_tree t

(** Calculate the simplest evolutionary tree for a set of four species as given
    by their helices. Raises: Invalid_argument if the given helices are not all
    the same length. *)
let find_simplest_tree (x1 : helix) (x2 : helix) (x3 : helix) (x4 : helix) :
    tree = 
    all_trees x1 x2 x3 x4 |> add_ancestor_labels_list |> simplest_tree |> fst |> unlabel_tree

(** Every triplet of adjacent nucleotides in a helix encodes one of 20 amino
    acids, or else indicates the END of the chain. The amino acids and END
    marker are named in the following enumeration. *)
type acid =
  | Ala
  | Arg
  | Asn
  | Asp
  | Cys
  | Glu
  | Gln
  | Gly
  | His
  | Ile
  | Leu
  | Lys
  | Met
  | Phe
  | Pro
  | Ser
  | Thr
  | Trp
  | Tyr
  | Val
  | END

let string_of_acid = function
  | Ala -> "Ala"
  | Arg -> "Arg"
  | Asn -> "Asn"
  | Asp -> "Asp"
  | Cys -> "Cys"
  | Glu -> "Glu"
  | Gln -> "Gln"
  | Gly -> "Gly"
  | His -> "His"
  | Ile -> "Ile"
  | Leu -> "Leu"
  | Lys -> "Lys"
  | Met -> "Met"
  | Phe -> "Phe"
  | Pro -> "Pro"
  | Ser -> "Ser"
  | Thr -> "Thr"
  | Trp -> "Trp"
  | Tyr -> "Tyr"
  | Val -> "Val"
  | END -> "END"

(** Outputs the encoded amino acid or END marker corresponding to a given
    nucleotide triplet. *)
let acid_of_triplet (n1 : nucleotide) (n2 : nucleotide) (n3 : nucleotide) : acid
    =
  match (n1, n2, n3) with
  | A, A, A -> Phe
  | A, A, G -> Phe
  | A, A, T -> Leu
  | A, A, C -> Leu
  | G, A, A -> Leu
  | G, A, G -> Leu
  | G, A, T -> Leu
  | G, A, C -> Leu
  | T, A, A -> Ile
  | T, A, G -> Ile
  | T, A, T -> Ile
  | T, A, C -> Met
  | C, A, A -> Val
  | C, A, G -> Val
  | C, A, T -> Val
  | C, A, C -> Val
  | A, G, A -> Ser
  | A, G, G -> Ser
  | A, G, T -> Ser
  | A, G, C -> Ser
  | G, G, A -> Pro
  | G, G, G -> Pro
  | G, G, T -> Pro
  | G, G, C -> Pro
  | T, G, A -> Thr
  | T, G, G -> Thr
  | T, G, T -> Thr
  | T, G, C -> Thr
  | C, G, A -> Ala
  | C, G, G -> Ala
  | C, G, T -> Ala
  | C, G, C -> Ala
  | A, T, A -> Tyr
  | A, T, G -> Tyr
  | A, T, T -> END
  | A, T, C -> END
  | G, T, A -> His
  | G, T, G -> His
  | G, T, T -> Gln
  | G, T, C -> Gln
  | T, T, A -> Asn
  | T, T, G -> Asn
  | T, T, T -> Lys
  | T, T, C -> Lys
  | C, T, A -> Asp
  | C, T, G -> Asp
  | C, T, T -> Glu
  | C, T, C -> Glu
  | A, C, A -> Cys
  | A, C, G -> Cys
  | A, C, T -> END
  | A, C, C -> Trp
  | G, C, A -> Arg
  | G, C, G -> Arg
  | G, C, T -> Arg
  | G, C, C -> Arg
  | T, C, A -> Ser
  | T, C, G -> Ser
  | T, C, T -> Arg
  | T, C, C -> Arg
  | C, C, A -> Gly
  | C, C, G -> Gly
  | C, C, T -> Gly
  | C, C, C -> Gly

(** Given a helix, decodes its first acid chain, according to the following
    steps.

    1. Start decoding by scanning the helix for the first occurrence of the
    triplet TAC, which encodes the amino acid Met (and is the only triplet that
    encodes Met). This signals the start of an acid chain with Met as its head
    element.

    2. After the first occurrence of TAC, decode the subsequent nucleotides in
    the helix in triplets and build a list of the results.

    3. Stop decoding when either (a) there are fewer than three nucleotides
    remaining in the list, or (b) one of the triplets encoding the END acid
    (ATT, ACT, or ATC) is encountered. If we encounter the END acid, the END
    marker should NOT be included in the output acid list.

    Any nucleotides that occur before the first occurrence of Met, including END
    markers, are disregarded. If there are no occurrences of Met, the empty list
    is returned. *)
let rec acids_of_helix (x : helix) : acid list = 
  match x with 
  |[] -> []
  |[a1] -> []
  |a1 :: [a2] -> []
  |a1 :: a2 :: a3 :: t1 -> if acid_of_triplet a1 a2 a3 = Met 
    then match t1 with
    |[] -> []
    |[b1] -> []
    |b1 :: [b2] -> []
    |b1 :: b2 :: b3 :: t2 -> if acid_of_triplet b1 b2 b3 != END then
      acid_of_triplet b1 b2 b3 :: acids_of_helix t2 else []
     else match x with 
     |h::t -> acids_of_helix t

(** Decodes all the acid chains in a given helix, according to the following
    rules. The next acid chain starts after the END triplet of the previous
    chain. The beginning marker of the next chain (TAC) may not be immediately
    after the end of the previous one; there could be any number of nucleotides
    in between, and not necessarily a multiple of three. The chains do not
    overlap; so if in the middle of one amino acid-encoding nucleotide sequence
    there is another TAC, it is ignored. *)
let rec all_acids_of_helix (x : helix) : acid list list =
  failwith "Unimplemented"

(**********************************************************************)
(******************** functions to display results ********************)
(**********************************************************************)

(* Computes an "ascii art" representation of an evolutionary tree to display on
   the console. *)
let rec string_of_tree (r : tree) : string =
  let spaces (n : int) : string = String.make n ' ' in
  let dashes (n : int) : string = String.make n '-' in

  let rec zip ((t1, w1) : string list * int) ((t2, w2) : string list * int) :
      string list * int =
    let rec aux (t1 : string list) (t2 : string list) : string list =
      begin
        match (t1, t2) with
        | [], _ -> List.map (fun s -> spaces w1 ^ s) t2
        | _, [] -> List.map (fun s -> s ^ spaces w2) t1
        | t :: ts, s :: ss -> (t ^ s) :: aux ts ss
      end
    in
    let h1 = w1 / 2 in
    let h2 = w2 / 2 in
    let line1 = spaces h1 ^ dashes (h1 + h2 + 2) ^ spaces h2 in
    let line2 = spaces h1 ^ "|" ^ spaces (h1 + h2) ^ "|" ^ spaces h2 in
    (line1 :: line2 :: aux t1 t2, w1 + w2 + 2)
  in

  let rec aux (t : tree) : string list * int =
    begin
      match t with
      | Leaf x ->
          let str1 = " " ^ string_of_ape x ^ " " in
          let str2 =
            if String.length str1 mod 2 = 0 then str1 else " " ^ str1
          in
          ([ str2 ], String.length str2)
      | Node (left, right) -> zip (aux left) (aux right)
    end
  in
  let strs, w = aux r in
  let strs2 = (spaces (w / 2) ^ " |" ^ spaces (w / 2)) :: strs in
  String.concat "\n" strs2 ^ "\n"

(* Creates an "ascii art" representation of all trees in a list of trees for
   display on the console. *)
let rec string_of_tree_list (ts : tree list) : string =
  begin
    match ts with
    | [] -> ""
    | hd :: tl -> string_of_tree hd ^ "\n\n" ^ string_of_tree_list tl
  end

(* Displays the number of mutations between human and each ape. *)
let display_mutations () =
  let mutations (ape : helix) : string =
    Printf.sprintf "Number of mutations between humans and %s: %d"
      (string_of_ape ape)
      (hamming_distance human ape)
  in
  print_endline (mutations gorilla);
  print_endline (mutations lar_gibbon);
  print_endline (mutations pileated_gibbon);
  print_endline (mutations siamang);
  print_endline (mutations white_cheeked_gibbon);
  print_endline (mutations orangutan);
  print_endline (mutations chimpanzee)

(* Displays the greater ape subtree in Figure A. *)
let display_greater_ape_tree () =
  print_endline "Greater apes (from Figure A)";
  print_endline "----------------------------";
  print_endline (string_of_tree (greater_apes ()))

(* Displays the lesser ape subtree in Figure A. *)
let display_lesser_ape_tree () =
  print_endline "Lesser apes (from Figure A)";
  print_endline "---------------------------";
  print_endline (string_of_tree (lesser_apes ()))

(* Displays Figure A. *)
let display_figure_A () =
  print_endline "Figure A";
  print_endline "--------";
  print_endline (string_of_tree (Node (greater_apes (), lesser_apes ())))

(* Displays all evolutionary trees of the greater apes. *)
let display_greater_ape_trees () =
  print_endline "Possible evolutionary trees (greater apes)";
  print_endline "------------------------------------------";
  print_endline (string_of_tree_list (all_greater_ape_trees ()))

(* Displays all evolutionary trees of the lesser apes. *)
let display_lesser_ape_trees () =
  print_endline "Possible evolutionary trees (lesser apes)";
  print_endline "------------------------------------------";
  print_endline (string_of_tree_list (all_lesser_ape_trees ()))

(** Displays the simplest evolutionary tree for the greater and lesser apes. *)
let display_simplest_trees () : unit =
  print_endline "Computed evolutionary tree for greater apes";
  print_endline "-------------------------------------------";
  print_endline (string_of_tree (simplest_greater_ape_tree ()));
  print_endline "";
  print_endline "Computed evolutionary tree for lesser apes";
  print_endline "------------------------------------------";
  print_endline (string_of_tree (simplest_lesser_ape_tree ()))

let hours_worked = -1
