/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.util;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.concurrent.LinkedBlockingQueue;

/**
 *@author Pavel Silin, Pavel_Silin@epam.com
 **/

public class BatchBuffer<T> {

    public static final SAMRecord FINISH_RECORD = new SAMRecord(null);
    private Collection<T> batch;
    private int bufferSize;
    private LinkedBlockingQueue<Collection<T>> batchQueue;

    public BatchBuffer(int bufferSize){
        this.bufferSize = bufferSize;
        batch = new ArrayList<T>(bufferSize);
        batchQueue = new LinkedBlockingQueue<Collection<T>>();
    }

    public BatchBuffer(int bufferSize, int batchCapacity){
        this.bufferSize = bufferSize;
        batch = new ArrayList<T>(bufferSize);
        batchQueue = new LinkedBlockingQueue<Collection<T>>(batchCapacity);
    }


    public void add(T element) throws InterruptedException {
        batch.add(element);

        if (batch.size() == bufferSize || element == FINISH_RECORD) {
            batchQueue.put(batch);
            batch = buildBatch();
        }
    }

    public Collection<T> take() throws InterruptedException {
        return batchQueue.take();
    }

    protected Collection<T> buildBatch() {
        return new ArrayList<T>(bufferSize);
    }
}
