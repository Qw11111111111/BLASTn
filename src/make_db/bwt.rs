use std::rc::Rc;
#[derive(Default, Debug)]
struct Bwt {
    _last_col: Vec<Rc<char>>,
    _null_rot: usize,
}

pub fn bwt(data: &str) -> Result<(), String> {
    let _bwt = gen_bwt(data);

    Ok(())
}

fn _save_bwt_idx() -> Result<(), String> {
    Ok(())
}

fn _get_bwt_idx() {}

fn gen_bwt(_data: &str) -> Bwt {
    Bwt::default()
}
