use std::{
    cmp::{Eq, Ordering, PartialEq},
    fmt::{Display, Write},
};

type Array = Vec<Vec<i32>>;
type PathArrray = Vec<Vec<Vec<Move>>>;
type GapArry = Vec<Vec<Gap>>;

#[derive(Debug, Clone)]
enum Move {
    Down,
    Right,
    Diag,
}

#[derive(Debug, Clone, Eq, PartialEq)]
enum Gap {
    Gap,
    NoGap,
}

pub struct Arguments {
    pub db: String,
    pub query: String,
    pub local: bool,
    pub verbose: bool,
}

#[derive(Default, Debug)]
struct ScoreMatrix {
    scores: Array,
    paths: PathArrray,
    gaps: GapArry,
}

impl ScoreMatrix {
    fn init(&mut self, seq1: &str, seq2: &str) {
        //           rows   cols
        self.scores = vec![vec![0; seq2.len() + 1]; seq1.len() + 1];
        self.paths = vec![vec![Vec::default(); seq2.len() + 1]; seq1.len() + 1];
        self.gaps = vec![vec![Gap::NoGap; seq2.len() + 1]; seq1.len() + 1];
    }

    fn fill(&mut self, seq1: &str, seq2: &str, scheme: &ScoringScheme, local: bool) {
        if !local {
            self.scores[0]
                .iter_mut()
                .zip(self.paths[0].iter_mut())
                .zip(self.gaps[0].iter_mut())
                .enumerate()
                .for_each(|(i, ((item, path), gap))| {
                    *item += i as i32 * scheme.gap_extension + scheme.gap_opening;
                    path.append(&mut vec![Move::Right]);
                    *gap = Gap::Gap;
                });
            self.scores
                .iter_mut()
                .zip(self.paths.iter_mut())
                .zip(self.gaps.iter_mut())
                .enumerate()
                .for_each(|(i, ((row, path_row), gap))| {
                    row[0] += i as i32 * scheme.gap_extension + scheme.gap_opening;
                    path_row[0].append(&mut vec![Move::Down]);
                    gap[0] = Gap::Gap;
                });
        }
        for i in 1..=seq1.len() {
            for j in 1..=seq2.len() {
                let diag_score = if seq1.chars().nth(i - 1) == seq2.chars().nth(j - 1) {
                    self.scores[i - 1][j - 1] + scheme.match_
                } else {
                    self.scores[i - 1][j - 1] + scheme.mismatch
                };
                let down_score = if self.gaps[i - 1][j] == Gap::Gap {
                    self.scores[i - 1][j] + scheme.gap_extension
                } else {
                    self.scores[i - 1][j] + scheme.gap_opening
                };
                let right_score = if self.gaps[i][j - 1] == Gap::Gap {
                    self.scores[i][j - 1] + scheme.gap_extension
                } else {
                    self.scores[i][j - 1] + scheme.gap_opening
                };
                let max_score = down_score.max(right_score).max(diag_score);
                //let p: Vec<Move> = self.paths[i - 1][j].clone();
                if max_score == down_score || max_score == right_score {
                    self.gaps[i][j] = Gap::Gap;
                }
                if local && max_score < 0 {
                    self.paths[i][j] = Vec::default();
                } else {
                    self.scores[i][j] = max_score;
                    if max_score == down_score {
                        self.paths[i][j].push(Move::Down);
                    }
                    if max_score == right_score {
                        self.paths[i][j].push(Move::Right);
                    }
                    if max_score == diag_score {
                        self.paths[i][j].push(Move::Diag);
                    }
                }
            }
        }
    }

    fn backtrace(&self, seq1: &str, seq2: &str, local: bool) {
        let start: Vec<(usize, usize)> = if local {
            argmax(&self.scores)
        } else {
            vec![(seq1.len(), seq2.len())]
        };
        for start_ in &start {
            let mut hit = Hit::default();
            if !get_next(&self.paths, seq1, seq2, *start_, &mut hit) {
                println!("could not generate alignement");
            }
        }
    }
}

impl Display for ScoreMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Scores: ")?;
        for row in &self.scores {
            for item in row {
                write!(f, "{} ", item)?;
            }
            writeln!(f)?;
        }
        Ok(())
        /*    writeln!(f, "paths: ")?;
            for row in &self.paths {
                for p in row {
                    write!(f, "{:?}", p)?;
                }
                writeln!(f)?;
            }
            Ok(())
        }*/
    }
}

#[derive(Default, Debug)]
struct ScoringScheme {
    gap_opening: i32,
    gap_extension: i32,
    mismatch: i32,
    match_: i32,
}

#[derive(Default, Debug)]
struct Hit {
    query: String,
    db: String,
    start_in_query: usize,
    start_in_db: usize,
}

impl Display for Hit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "\nseq1: {}",
            self.query.chars().rev().collect::<String>()
        )?;
        write!(f, "\n      ")?;
        for (s1, s2) in self.query.chars().rev().zip(self.db.chars().rev()) {
            if s1 == s2 {
                write!(f, "|")?;
            } else {
                write!(f, " ")?;
            }
        }
        write!(f, "\nseq2: {}\n", self.db.chars().rev().collect::<String>())?;
        writeln!(
            f,
            "start in seq1: {}\nstart in seq2: {}\n",
            self.start_in_query, self.start_in_db
        )?;
        Ok(())
    }
}

pub fn needleman_wunsch(args: &Arguments) {
    let scheme = ScoringScheme {
        gap_opening: -8,
        gap_extension: -6,
        mismatch: -4,
        match_: 5,
    };
    let mut mat = ScoreMatrix::default();
    mat.init(&args.query, &args.db);
    mat.fill(&args.query, &args.db, &scheme, args.local);
    if args.verbose {
        println!("{}", mat);
    }
    mat.backtrace(&args.query, &args.db, args.local);
}

fn get_next(
    paths: &PathArrray,
    seq1: &str,
    seq2: &str,
    current: (usize, usize),
    hit: &mut Hit,
) -> bool {
    if current == (0, 0) {
        //println!("\nHit: {}\n", hit);
        return true;
    }
    if paths[current.0][current.1].is_empty() {
        //println!("\nHit: {}\n", hit);
        return true;
    }
    for p in &paths[current.0][current.1] {
        hit.start_in_query = current.0.max(1) - 1;
        hit.start_in_db = current.1.max(1) - 1;
        let next: (usize, usize) = match p {
            Move::Down => {
                hit.query
                    .write_char(seq1.chars().nth(current.0 - 1).unwrap_or_else(|| {
                        panic!("could not index into seq1 at {}", current.0 - 1)
                    }))
                    .expect("could not append char");
                hit.db.write_char('-').expect("could not append char");
                (current.0 - 1, current.1)
            }
            Move::Right => {
                hit.query.write_char('-').expect("could not append char");
                hit.db
                    .write_char(seq2.chars().nth(current.1 - 1).unwrap_or_else(|| {
                        panic!("could not index into seq2 at {}", current.1 - 1)
                    }))
                    .expect("could not append char");
                (current.0, current.1 - 1)
            }
            Move::Diag => {
                hit.query
                    .write_char(seq1.chars().nth(current.0 - 1).unwrap_or_else(|| {
                        panic!("could not index into seq1 at {}", current.0 - 1)
                    }))
                    .expect("could not append char");
                hit.db
                    .write_char(seq2.chars().nth(current.1 - 1).unwrap_or_else(|| {
                        panic!("could not index into seq2 at {}", current.1 - 1)
                    }))
                    .expect("could not append char");
                (current.0 - 1, current.1 - 1)
            }
        };
        if get_next(paths, seq1, seq2, next, hit) {
            return true;
        }
        hit.query.pop();
        hit.db.pop();
    }
    false
}

fn argmax(mat: &Array) -> Vec<(usize, usize)> {
    let mut idc = Vec::default();
    let mut max = i32::MIN;
    for i in 0..mat.len() {
        for j in 0..mat[0].len() {
            match mat[i][j].cmp(&max) {
                Ordering::Greater => {
                    max = mat[i][j];
                    idc = vec![(i, j)];
                }
                Ordering::Equal => idc.push((i, j)),
                _ => (),
            }
        }
    }
    idc
}
