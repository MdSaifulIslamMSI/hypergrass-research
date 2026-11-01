package hypergrass;

/** Diagnostic: how do single-byte values survive the range coder? */
public final class RangeCoderDiag {
    public static void main(String[] args) {
        int ok = 0;
        for (int v = 0; v < 256; v++) {
            byte[] data = {(byte) v};
            byte[] comp = RangeCoder.encode(data, 1);
            byte[] back = RangeCoder.decode(comp, 1);
            if (back[0] == (byte) v) {
                ok++;
            } else if (v < 8) {
                System.out.printf("%d -> %d (comp %02x %02x %02x %02x %02x)%n",
                        v, back[0] & 0xFF, comp[0] & 0xFF, comp[1] & 0xFF, comp[2] & 0xFF,
                        comp[3] & 0xFF, comp[4] & 0xFF);
            }
        }
        System.out.println("single-byte round-trips: " + ok + "/256");
    }
}
