pub fn arange(start: f64, stop: f64, step: usize) -> impl Iterator<Item = f64> {
    let step = step - 1;
    let mut count = start;
    let ediv = (stop - start) / step as f64;
    let mut done = false;

    std::iter::from_fn(move || {
        if count <= stop {
            let out = Some(count);
            count += ediv;

            out
        } else {
            None
        }
    })
}
